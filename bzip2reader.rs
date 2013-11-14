use super::checksum::StrangeCRC;
use super::consts;

use std::io::{ Reader};
use std::vec;
use std::char;
use std::iter::range_inclusive;

#[allow(missing_doc)]
pub struct Bzip2Reader<R> {
    /// Index of last char in the block, so the block size == last + 1
    priv last: i32,
    /// Index in zptr[] of original string before sorting
    orig_ptr: i32,

    /// always in range 0..9,
    /// Current block size = 100,000 * this number
    block_size100k: i32,

    block_randomised: bool,

    bytes_out: i32,
    bs_buff: i32,
    bs_live: i32,
    m_crc: StrangeCRC,

    in_use: ~[bool, ..256],
    n_inuse: i32,

    seq_to_unseq: ~[u8, ..256],
    unseq_to_seq: ~[u8, ..256],

    selector: ~[u8, ..consts::MAXIMUM_SELECTORS],
    selector_mtf: ~[u8, ..consts::MAXIMUM_SELECTORS],

    tt: ~[i32],
    ll8: ~[u8],

    /// freq table collected to save a pass over the data
    /// during decompression.
    unzftab: ~[i32, ..256],

    limit: ~[~[i32]],
    base_array: ~[~[i32]],
    perm: ~[~[i32]],
    min_lens: ~[i32, ..consts::GROUP_COUNT],

    base_reader: R,
    stream_end: bool,

    current_char: i32,
    current_state: State,

    stored_block_crc: i32,
    stored_combined_crc: i32,
    computed_block_crc: i32,
    computed_combined_crc: u32,

    count: i32,
    chprev: i32,
    ch2: i32,
    tPos: i32,
    rNToGo: i32,
    rTPos: i32,
    i2: i32,
    j2: i32,
    z: u8,
    is_owner: bool
}

impl<R: Reader> Reader for Bzip2Reader<R> {
    fn eof(&mut self) -> bool {
        self.stream_end || self.base_reader.eof()
    }

    fn read(&mut self, buf: &mut[u8]) -> Option<uint> {
        let mut ret = buf.len();

        for i in range(0, buf.len()) {
            let rb_res = self.read_b();
            match rb_res {
                Some(byte) => { buf[i] = byte; }
                None => { ret =  i; }
            };
        }

        match ret {
            0 => None,
            x => Some(x)
        }
    }
}

impl<R: Reader> Bzip2Reader<R> {
    /// Create a new stream reader
    /// # Arguments
    /// * `reader` - underlying reader to decompress from
    pub fn new(reader: R) -> Bzip2Reader<R> {
        let mut br = Bzip2Reader {
            last: 0,
            orig_ptr: 0,
            block_size100k: 0,
            block_randomised: false,
            bytes_out: 0,
            bs_buff: 0,
            bs_live: 0,
            m_crc: StrangeCRC::new(),
            in_use: ~([false, ..256]),
            n_inuse: 0,
            seq_to_unseq: ~([0, ..256]),
            unseq_to_seq: ~([0, ..256]),
            selector: ~([0, ..consts::MAXIMUM_SELECTORS]),
            selector_mtf: ~([0, ..consts::MAXIMUM_SELECTORS]),
            tt: ~[],
            ll8: ~[],
            unzftab: ~([0, ..256]),
            limit: vec::from_fn(consts::GROUP_COUNT as uint,
                |_| vec::from_elem(consts::MAXIMUM_ALPHA_SIZE as uint, 0i32)),
            base_array: vec::from_fn(consts::GROUP_COUNT as uint,
                |_| vec::from_elem(consts::MAXIMUM_ALPHA_SIZE as uint, 0i32)),
            perm: vec::from_fn(consts::GROUP_COUNT as uint,
                |_| vec::from_elem(consts::MAXIMUM_ALPHA_SIZE as uint, 0i32)),
            min_lens: ~([0, ..consts::GROUP_COUNT]),
            base_reader: reader,
            stream_end: false,
            current_char: 0,
            current_state: START_BLOCK_STATE,
            stored_block_crc: 0,
            stored_combined_crc: 0,
            computed_block_crc: 0,
            computed_combined_crc: 0,
            count: 0,
            chprev: 0,
            ch2: 0,
            tPos: 0,
            rNToGo: 0,
            rTPos: 0,
            i2: 0,
            j2: 0,
            z: 0,
            is_owner: false
        };

        br.bs_set_stream();
        br.initialize();
        br.init_block();
        br.setup_block();
        br
    }

    #[inline(always)]
    fn read_b(&mut self) -> Option<u8> {
        if self.stream_end {
            return None;
        }

        let ret: i32 = self.current_char;
        match self.current_state {
            RAND_PART_B_STATE => self.setup_randB(),
            RAND_PART_C_STATE => self.setup_randC(),
            NO_RAND_PART_B_STATE => self.setup_norandB(),
            NO_RAND_PART_C_STATE => self.setup_norandC(),
            START_BLOCK_STATE | NO_RAND_PART_A_STATE | RAND_PART_A_STATE => {}
        }

        Some(ret as u8)
    }

    fn make_maps(&mut self) {
        self.n_inuse = 0;

        for i in range(0, 256) {
            if self.in_use[i] {
                self.seq_to_unseq[self.n_inuse as uint] = i as u8;
                self.unseq_to_seq[i] = self.n_inuse as u8;
                self.n_inuse += 1;
            }
        }
    }

    fn initialize(&mut self) {
        let magic1: char = self.bs_get_uchar();
        let magic2: char = self.bs_get_uchar();
        let magic3: char = self.bs_get_uchar();
        let magic4: char = self.bs_get_uchar();

        if magic1 != 'B' || magic2 != 'Z' || magic3 != 'h' || magic4 < '1' || magic4 > '9' {
            self.stream_end = true;
            return;
        }

        self.set_decompress_structure_sizes((magic4 as u8 - '0' as u8) as i32);
        self.computed_combined_crc = 0;
    }

    fn init_block(&mut self) {
        let magic1: char = self.bs_get_uchar();
        let magic2: char = self.bs_get_uchar();
        let magic3: char = self.bs_get_uchar();
        let magic4: char = self.bs_get_uchar();
        let magic5: char = self.bs_get_uchar();
        let magic6: char = self.bs_get_uchar();

        if magic1 == '\x17' && magic2 == '\x72' && magic3 == '\x45' && magic4 == '\x38'
            && magic5 == '\x50' && magic6 == '\x90' {
            self.complete();
            return;
        }

        if magic1 != '\x31' || magic2 != '\x41' || magic3 != '\x59' || magic4 != '\x26'
            || magic5 != '\x53' || magic6 != '\x59' {
            fail!("Bad Bzip2 block header in `init_block`");
//             self.stream_end = true;
        }

        self.stored_block_crc = self.bs_get_int32();
        self.block_randomised = self.bs_r(1) == 1;

        self.get_and_move_to_front_decode();
        self.m_crc.reset();
        self.current_state = START_BLOCK_STATE;
    }

    fn end_block(&mut self) {
        self.computed_block_crc = self.m_crc.value();

        // Bad CRC is fatal
        if self.stored_block_crc != self.computed_block_crc {
            fail!("Fatal CRC error in `end_block`");
        }

        // 1528150659
        self.computed_combined_crc = ((self.computed_combined_crc << 1) & 0xFF_FF_FF_FF) |
            (self.computed_combined_crc >> 31);
        self.computed_combined_crc = self.computed_combined_crc ^ self.computed_block_crc as u32;
    }

    fn complete(&mut self) {
        self.stored_combined_crc = self.bs_get_int32();
        if self.stored_combined_crc as u32 != self.computed_combined_crc {
            fail!("BZip2 input stream CRC error");
        }

        self.stream_end = true;
    }

    fn bs_set_stream(&mut self) {
//         self.base_reader = self.base_reader;
        self.bs_live = 0;
        self.bs_buff = 0;
    }

    fn fill_buffer(&mut self) {
        let thech: i32 = self.base_reader.read_byte()
            .expect("fill_buffer() failed to read byte") as i32;

        self.bs_buff = (self.bs_buff << 8) | (thech & 0xFF);
        self.bs_live += 8;
    }

    fn bs_r(&mut self, n: i32) -> i32 {
        while self.bs_live < n {
            self.fill_buffer();
        }

        let v: i32 = (self.bs_buff >> (self.bs_live - n)) & ((1 << n) - 1);
        self.bs_live -= n;
        v
    }

    fn bs_get_uchar(&mut self) -> char {
        char::from_u32(self.bs_r(8) as u32).expect("bs_r(8) returned invalid u32")
    }

    fn bs_get_intVS(&mut self, num_bits: i32) -> i32 {
        self.bs_r(num_bits)
    }

    fn bs_get_int32(&mut self) -> i32 {
        let mut res = self.bs_r(8);
        res = (res << 8) | self.bs_r(8);
        res = (res << 8) | self.bs_r(8);
        res = (res << 8) | self.bs_r(8);

        res
    }

    fn recv_decoding_tables(&mut self) {
        let mut len: ~[~[char]] = vec::from_fn(consts::GROUP_COUNT as uint,
            |_| vec::from_elem(consts::MAXIMUM_ALPHA_SIZE as uint, 0u8 as char));

        let mut in_use16: ~[bool] = vec::from_elem(16, false);

        // Receive the mapping table
        for i in range(0, 16) {
            in_use16[i] = self.bs_r(1) == 1;
        }

        for i in range(0, 16) {
            if in_use16[i] {
                for j in range(0, 16) {
                    self.in_use[i * 16 + j] = self.bs_r(1) == 1;
                }
            } else {
                for j in range(0, 16) {
                    self.in_use[i * 16 + j] = false;
                }
            }
        }

        self.make_maps();
        let alpha_size: i32 = self.n_inuse + 2;

        // Now the selectors
        let n_groups: i32 = self.bs_r(3);
        let n_selectors: i32 = self.bs_r(15);

        for i in range(0, n_selectors) {
            let mut j: i32 = 0;
            while self.bs_r(1) == 1 {
                j += 1;
            }
            self.selector_mtf[i] = j as u8;
        }

        // Undo the MTF values for the selectors
        let mut pos: ~[u8] = vec::from_elem(consts::GROUP_COUNT as uint, 0u8);
        for v in range(0, n_groups) {
            pos[v] = v as u8;
        }

        for i in range(0, n_selectors) {
            let mut v: i32 = self.selector_mtf[i] as i32;
            let tmp = pos[v];
            while v > 0 {
                pos[v] = pos[v as uint - 1];
                v -= 1;
            }
            pos[0] = tmp;
            self.selector[i] = tmp;
        }

        // Now the coding tables
        for t in range(0, n_groups) {
            let mut curr: i32 = self.bs_r(5);
            for i in range(0, alpha_size) {
                while self.bs_r(1) == 1 {
                    if self.bs_r(1) == 0 {
                        curr += 1;
                    } else {
                        curr -= 1;
                    }
                }
                len[t][i] = char::from_u32(curr as u32).expect("error in recv_decoding_tables");
            }
        }

        // Create the Huffman decoding tables
        for t in range(0, n_groups) {
            let mut min_len: i32 = 32;
            let mut max_len: i32 = 0;
            for i in range(0, alpha_size) {
                max_len = ::std::num::max(max_len, len[t][i] as i32);
                min_len = ::std::num::min(min_len, len[t][i] as i32);
            }
            hb_create_decode_table(self.limit[t], self.base_array[t],
                self.perm[t], len[t], min_len, max_len, alpha_size);
            self.min_lens[t] = min_len;
        }
    }

    fn get_and_move_to_front_decode(&mut self) {
        let mut yy: [u8, ..256] = [0u8, ..256];
        let mut next_sym: i32;

        let limit_last: i32 = consts::BASE_BLOCK_SIZE * self.block_size100k;
        self.orig_ptr = self.bs_get_intVS(24);

        self.recv_decoding_tables();
        let EOB: i32 = self.n_inuse + 1;
        let mut group_no: i32 = -1;
        let mut group_pos: i32 = 0;

        // Setting up unzftables for perf reasons
        for i in range_inclusive(0, 255) {
            self.unzftab[i] = 0;
        }

        for i in range_inclusive(0, 255) {
            yy[i] = i as u8;
        }

        self.last = -1;

        if group_pos == 0 {
            group_no += 1;
            group_pos = consts::GROUP_SIZE;
        }

        group_pos -= 1;
        let mut zt: i32 = self.selector[group_no] as i32;
        let mut zn: i32 = self.min_lens[zt];
        let mut zvec: i32 = self.bs_r(zn);
        let mut zj: i32;

        while zvec > self.limit[zt][zn] {
            if zn > 20 { // teh longest code
                fail!("Bzip2 data error in `get_and_move_to_front_decode`");
            }
            zn += 1;

            while self.bs_live < 1 {
                self.fill_buffer();
            }

            zj = (self.bs_buff >> (self.bs_live -1)) & 1;
            self.bs_live -= 1;
            zvec = (zvec << 1) | zj;
        }

        if zvec - self.base_array[zt][zn] < 0 ||
            zvec - self.base_array[zt][zn] >= consts::MAXIMUM_ALPHA_SIZE {
            fail!("Bzip2 data error in `get_and_move_to_front_decode`");
        }

        next_sym = self.perm[zt][zvec - self.base_array[zt][zn]];

        loop {
            if next_sym == EOB {
                break;
            }

            if next_sym == consts::RUN_A || next_sym == consts::RUN_B {
                let mut s: i32 = -1;
                let mut n: i32 = 1;

                do_while!({
                    if next_sym == consts::RUN_A {
                        s += (0 + 1) * n;
                    } else if next_sym == consts::RUN_B {
                        s += (1 + 1) * n;
                    }

                    n <<= 1;

                    if group_pos == 0 {
                        group_no += 1;
                        group_pos = consts::GROUP_SIZE;
                    }

                    group_pos -= 1;

                    zt = self.selector[group_no] as i32;
                    zn = self.min_lens[zt];
                    zvec = self.bs_r(zn);

                    while zvec > self.limit[zt][zn] {
                        zn += 1;
                        while self.bs_live < 1 {
                            self.fill_buffer();
                        }
                        zj = (self.bs_buff >> (self.bs_live - 1)) & 1;
                        self.bs_live -= 1;
                        zvec = (zvec << 1) | zj;
                    }
                    next_sym = self.perm[zt][zvec - self.base_array[zt][zn]];
                }, (next_sym == consts::RUN_A || next_sym == consts::RUN_B));

                s += 1;
                let ch: u8 = self.seq_to_unseq[yy[0]];
                self.unzftab[ch as uint] += s;

                while s > 0 {
                    self.last += 1;
                    self.ll8[self.last as uint] = ch;
                    s -= 1;
                }

                if self.last >= limit_last {
                    fail!("BZip2 input stream block overrun");
                }
            } else {
                self.last += 1;
                if self.last >= limit_last {
                    fail!("BZip2 input stream block overrun");
                }

                let tmp: u8 = yy[next_sym - 1];
                self.unzftab[self.seq_to_unseq[tmp as uint]] += 1;
                self.ll8[self.last as uint] = self.seq_to_unseq[tmp as uint];

                {
                    let mut j: i32 = next_sym - 1;
                    while j > 0 {
                        yy[j as uint] = yy[j as int - 1];
                        j -= 1;
                    }
                }

                yy[0] = tmp;

                if group_pos == 0 {
                    group_no += 1;
                    group_pos = consts::GROUP_SIZE;
                }

                group_pos -= 1;

                zt = self.selector[group_no] as i32;
                zn = self.min_lens[zt];
                zvec = self.bs_r(zn);
                while zvec > self.limit[zt][zn] {
                    zn += 1;
                    while self.bs_live < 1 {
                        self.fill_buffer();
                    }
                    zj = (self.bs_buff >> (self.bs_live - 1)) & 1;
                    self.bs_live -= 1;
                    zvec = (zvec << 1) | zj;
                }
                next_sym = self.perm[zt][zvec - self.base_array[zt][zn]];
                continue;
            }
        }
    }

    fn setup_block(&mut self) {
        let mut cftab: ~[i32] = vec::from_elem(257, 0i32);

        cftab[0] = 0;
        for (i, &x) in self.unzftab.iter().enumerate() {
            cftab[i + 1] = x;
        }

        for i in range_inclusive(1, 256) {
            cftab[i] += cftab[i - 1];
        }

        for i in range_inclusive(0, self.last) {
            let ch: u8 = self.ll8[i];
            self.tt[cftab[ch]] = i;
            cftab[ch] += 1;
        }

        let _ = cftab;

        self.tPos = self.tt[self.orig_ptr];

        self.count = 0;
        self.i2 = 0;
        self.ch2 = 256; // Not char and not EOF

        if self.block_randomised {
            self.rNToGo = 0;
            self.rTPos = 0;
            self.setup_randA();
        } else {
            self.setup_norandA();
        }
    }

    fn setup_randA(&mut self) {
        if self.i2 <= self.last {
            self.chprev = self.ch2;
            self.ch2 = self.ll8[self.tPos] as i32;
            self.tPos = self.tt[self.tPos];
            if self.rNToGo == 0 {
                self.rNToGo = consts::RANDOM_NUMBERS[self.rTPos];
                self.rTPos += 1;
                if self.rTPos == 512 {
                    self.rTPos = 0;
                }
            }
            self.rNToGo -= 1;
            self.ch2 ^= if self.rNToGo == 1 { 1 } else { 0 };
            self.i2 += 1;

            self.current_char = self.ch2;
            self.current_state = RAND_PART_B_STATE;
            self.m_crc.update_val(self.ch2);
        } else {
            self.end_block();
            self.init_block();
            self.setup_block();
        }
    }

    fn setup_norandA(&mut self) {
        if self.i2 <= self.last {
            self.chprev = self.ch2;
            self.ch2 = self.ll8[self.tPos] as i32;
            self.tPos = self.tt[self.tPos];
            self.i2 += 1;

            self.current_char = self.ch2;
            self.current_state = NO_RAND_PART_B_STATE;
            self.m_crc.update_val(self.ch2);
        } else {
            self.end_block();
            self.init_block();
            self.setup_block();
        }
    }

    fn setup_randB(&mut self) {
        if self.ch2 != self.chprev {
            self.current_state = RAND_PART_A_STATE;
            self.count = 1;
            self.setup_randA();
        } else {
            self.count += 1;
            if self.count >= 4 {
                self.z = self.ll8[self.tPos];
                self.tPos = self.tt[self.tPos];
                if self.rNToGo == 0 {
                    self.rNToGo = consts::RANDOM_NUMBERS[self.rTPos];
                    self.rTPos += 1;
                    if self.rTPos == 512 {
                        self.rTPos = 0;
                    }
                }
                self.rNToGo -= 1;
                self.z ^= if self.rNToGo == 1 { 1 } else { 0 } as u8;
                self.j2 = 0;
                self.current_state = RAND_PART_C_STATE;
                self.setup_randC();
            } else {
                self.current_state = RAND_PART_A_STATE;
                self.setup_randA();
            }
        }
    }

    fn setup_randC(&mut self) {
        if self.j2 < self.z as i32 {
            self.current_char = self.ch2;
            self.m_crc.update_val(self.ch2);
            self.j2 += 1;
        } else {
            self.current_state = RAND_PART_A_STATE;
            self.i2 += 1;
            self.count = 0;
            self.setup_randA();
        }
    }

    fn setup_norandB(&mut self) {
        if self.ch2 != self.chprev {
            self.current_state = NO_RAND_PART_A_STATE;
            self.count = 1;
            self.setup_norandA();
        } else {
            self.count += 1;
            if self.count >= 4 {
                self.z = self.ll8[self.tPos];
                self.tPos = self.tt[self.tPos];
                self.current_state = NO_RAND_PART_C_STATE;
                self.j2 = 0;
                self.setup_norandC();
            } else {
                self.current_state = NO_RAND_PART_A_STATE;
                self.setup_norandA();
            }
        }
    }

    fn setup_norandC(&mut self) {
        if self.j2 < self.z as i32 {
            self.current_char = self.ch2;
            self.m_crc.update_val(self.ch2);
            self.j2 += 1;
        } else {
            self.current_state = NO_RAND_PART_A_STATE;
            self.i2 += 1;
            self.count = 0;
            self.setup_norandA();
        }
    }

    fn set_decompress_structure_sizes(&mut self, new_size100k: i32) {
        if !(0 <= new_size100k && new_size100k <= 9 &&
            0 <= self.block_size100k && self.block_size100k <= 9) {
            fail!("Invalid blocksize");
        }

        self.block_size100k = new_size100k;
        if new_size100k == 0 {
            return;
        }

        let n: i32 = consts::BASE_BLOCK_SIZE * new_size100k;
        self.ll8 = vec::from_elem(n as uint, 0u8);
        self.tt = vec::from_elem(n as uint, 0i32);
    }
}

fn hb_create_decode_table(limit: &mut [i32],
    base_array: &mut [i32],
    perm: &mut [i32],
    length: &[char],
    min_len: i32,
    max_len: i32,
    alpha_size: i32) {

    let mut pp: i32 = 0;

    for i in range_inclusive(min_len, max_len) {
        for j in range(0, alpha_size) {
            if length[j] as i32 == i {
                perm[pp] = j;
                pp += 1;
            }
        }
    }

    for i in range(0, consts::MAXIMUM_CODE_LENGTH) {
        base_array[i] = 0;
    }

    for i in range(0, alpha_size) {
        base_array[length[i] as uint + 1] += 1;
    }

    for i in range(1, consts::MAXIMUM_CODE_LENGTH) {
        base_array[i] += base_array[i - 1];
    }

    for i in range(0, consts::MAXIMUM_CODE_LENGTH) {
        limit[i] = 0;
    }

    let mut vec: i32 = 0;

    for i in range_inclusive(min_len, max_len) {
        vec += base_array[i + 1] - base_array[i];
        limit[i] = vec - 1;
        vec <<= 1;
    }

    for i in range_inclusive(min_len + 1, max_len) {
            base_array[i] = ((limit[i - 1] + 1) << 1) - base_array[i];
    }
}

enum State {
    START_BLOCK_STATE = 1,
    RAND_PART_A_STATE,
    RAND_PART_B_STATE,
    RAND_PART_C_STATE,
    NO_RAND_PART_A_STATE,
    NO_RAND_PART_B_STATE,
    NO_RAND_PART_C_STATE,
}

#[cfg(test)]
mod test {
    use super::Bzip2Reader;
    use std::io::{File};
    use std::vec;

//     #[test]
    fn test_correct() {
        let orig_path = Path::new("LICENSE.txt");
        let comp_path = Path::new("LICENSE.txt.bz2");

        let mut orig_f = File::open(&orig_path);
        let comp_f = File::open(&comp_path);

        let mut bzr = Bzip2Reader::new(comp_f);

        let orig_data = orig_f.read_to_end();
        let decompressed_data = bzr.read_to_end();

        if orig_data != decompressed_data {
            let os = ::std::str::from_utf8(orig_data);
            let ds = ::std::str::from_utf8(decompressed_data);
            println!("Original data: \n[{:s}]\n Decompressed data: \n[{:s}]", os, ds);
            fail!("Original data != Decompressed data");
        }
    }

    #[test]
    fn other_test() {
        let orig_path = Path::new("LICENSE.txt");
        let comp_path = Path::new("LICENSE.txt.bz2");

        let mut orig_f = File::open(&orig_path);
        let comp_f = File::open(&comp_path);

        let mut bzr = Bzip2Reader::new(comp_f);

        let orig_data = orig_f.read_to_end();

        let mut buffer = [0u8, ..1024 * 32];
        let mut decompressed_data = ~[];
        while !bzr.eof() {
            let n_read = bzr.read(buffer);
            match n_read {
                Some(n) => {
                    decompressed_data = vec::append(decompressed_data, buffer.slice_to(n))
                }
                None => { fail!("AAAAAH!"); }
            }
        }

//         let decompressed_data = bzr.read_to_end();

        if orig_data != decompressed_data {
            let os = ::std::str::from_utf8(orig_data);
            let ds = ::std::str::from_utf8(decompressed_data);
            println!("Original data: \n[{:s}]\n Decompressed data: \n[{:s}]", os, ds);
            fail!("Original data != Decompressed data");
        }
    }
}