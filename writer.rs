use checksum::StrangeCRC;
use consts;

use std::io::Writer;
use std::iter::range_inclusive;
use std::vec;
use std::char::from_u32;

pub struct Bzip2Writer<W> {
    last: i32,

    orig_ptr: i32,

    block_size100k: i32,

    block_randomised: bool,

    bytes_out: i32,
    bs_buff: i32,
    bs_live: i32,
    m_crc: StrangeCRC,

    in_use: ~[bool, ..256],
    n_inuse: i32,

    seq_to_unseq: ~[char, ..256],
    unseq_to_seq: ~[char, ..256],

    selector: ~[u8, ..consts::MAXIMUM_SELECTORS],
    selector_mtf: ~[u8, ..consts::MAXIMUM_SELECTORS],

    block: ~[u8],
    quadrant: ~[i32],
    zptr: ~[i32],
    szptr: ~[i16],
    ftab: ~[i32],

    n_mtf: i32,

    mt_freq: ~[i32, ..consts::MAXIMUM_ALPHA_SIZE],

    work_factor: i32,
    work_done: i32,
    work_limit: i32,
    first_attempt: bool,
    n_blocks_randomised: i32,

    current_char: i32,
    run_length: i32,
    blockcrc: u32,
    combinedcrc: u32,
    allowable_block_size: i32,

    base_stream: W,
}

impl<W: Writer> Drop for Bzip2Writer<W> {
    fn drop(&mut self) {
        if self.run_length > 0 {
            self.write_run();
        }

        self.current_char = -1;
        self.end_block();
        self.end_compression();
        self.flush();
    }
}

impl<W: Writer> Bzip2Writer<W> {
    pub fn new(stream: W, mut block_size: i32) -> Bzip2Writer<W> {
        if block_size < 1 {
            block_size = 1;
        }

        if block_size > 9 {
            block_size = 9;
        }

        let mut bw = Bzip2Writer {
            last: 0,
            orig_ptr: 0,
            block_size100k: block_size,
            block_randomised: false,
            bytes_out: 0,
            bs_buff: 0,
            bs_live: 0,
            m_crc: StrangeCRC::new(),
            in_use: ~([false, ..256]),
            n_inuse: 0,
            seq_to_unseq: ~([0u8 as char, ..256]),
            unseq_to_seq: ~([0u8 as char, ..256]),
            selector: ~([0, ..consts::MAXIMUM_SELECTORS]),
            selector_mtf: ~([0, ..consts::MAXIMUM_SELECTORS]),
            block: ~[],
            quadrant: ~[],
            zptr: ~[],
            szptr: ~[],
            ftab: ~[],
            n_mtf: 0,
            mt_freq: ~([0i32, ..consts::MAXIMUM_ALPHA_SIZE]),
            work_factor: 50,
            work_done: 0,
            work_limit: 0,
            first_attempt: false,
            n_blocks_randomised: 0,
            current_char: -1,
            run_length: 0,
            blockcrc: 0,
            combinedcrc: 0,
            allowable_block_size: 0,
            base_stream: stream
        };

        bw.bs_set_stream();
        bw.allocate_compression_structures();
        bw.initialize();
//         bw.init_block();

        bw
    }

    fn make_maps(&mut self) {
        self.n_inuse = 0;
        for i in range(0, 256) {
            if self.in_use[i] {
                self.seq_to_unseq[self.n_inuse] = ::std::char::from_u32(i as u32)
                    .expect("Error making char in `make_maps`");
                self.unseq_to_seq[i as int] = ::std::char::from_u32(self.n_inuse as u32)
                    .expect("Error making char in `make_maps`");
                self.n_inuse += 1;
            }
        }
    }

    fn write_run(&mut self) {
        if self.last < self.allowable_block_size {
            self.in_use[self.current_char] = true;
            for i in range(0, self.run_length) {
                self.m_crc.update_val(self.current_char);
            }

            match self.run_length {
                1 => {
                    self.last += 1;
                    self.block[self.last + 1] = self.current_char as u8;
                }
                2 => {
                    self.last += 1;
                    self.block[self.last + 1] = self.current_char as u8;
                    self.last += 1;
                    self.block[self.last + 1] = self.current_char as u8;
                }
                3 => {
                    self.last += 1;
                    self.block[self.last + 1] = self.current_char as u8;
                    self.last += 1;
                    self.block[self.last + 1] = self.current_char as u8;
                    self.last += 1;
                    self.block[self.last + 1] = self.current_char as u8;
                }
                _ => {
                    self.in_use[self.run_length - 4];
                    self.last += 1;
                    self.block[self.last + 1] = self.current_char as u8;
                    self.last += 1;
                    self.block[self.last + 1] = self.current_char as u8;
                    self.last += 1;
                    self.block[self.last + 1] = self.current_char as u8;
                    self.last += 1;
                    self.block[self.last + 1] = self.current_char as u8;
                    self.last += 1;
                    self.block[self.last + 1] = self.run_length as u8 - 4u8;
                }
            }
        } else {
            self.end_block();
            self.init_block();
        }
    }

    /// Get number of bytes written to the output
    pub fn bytes_written(&self) -> uint {
        self.bytes_out as uint
    }

    fn initialize(&mut self) {
        self.bytes_out = 0;
        self.n_blocks_randomised = 0;

        self.bs_put_uchar('B' as i32);
        self.bs_put_uchar('Z' as i32);
        self.bs_put_uchar('h' as i32);
        self.bs_put_uchar('0' as i32 + self.block_size100k as i32);

        self.combinedcrc = 0;
    }

    fn init_block(&mut self) {
        self.m_crc.reset();
        self.last = -1;

        for i in range(0, 256) {
            self.in_use[i] = false;
        }

        self.allowable_block_size = consts::BASE_BLOCK_SIZE * self.block_size100k - 20;
    }

    fn end_block(&mut self) {
        if self.last < 0 {
            return;
        }

        self.blockcrc = self.m_crc.value() as u32;
        self.combinedcrc = (self.combinedcrc << 1) | (self.combinedcrc >> 31);
        self.combinedcrc ^= self.blockcrc;

        // Sort the block and establish the position of the original string
        self.do_reversible_transformation();

        // 6 byte block header for damaged file recovery
        self.bs_put_uchar(0x31);
        self.bs_put_uchar(0x41);
        self.bs_put_uchar(0x59);
        self.bs_put_uchar(0x26);
        self.bs_put_uchar(0x53);
        self.bs_put_uchar(0x59);

        // Now the block's CRC
        self.bs_put_int(self.blockcrc as i32);

        // now a single bit indicating randomisation
        if self.block_randomised {
            self.bs_w(1, 1);
            self.n_blocks_randomised += 1;
        } else {
            self.bs_w(1, 0);
        }

        // Finally, block contents proper
        self.move_to_front_code_and_send();
    }

    fn end_compression(&mut self) {
        // 48-bit magic number sqrt(pi)
        self.bs_put_uchar(0x17);
        self.bs_put_uchar(0x72);
        self.bs_put_uchar(0x45);
        self.bs_put_uchar(0x38);
        self.bs_put_uchar(0x50);
        self.bs_put_uchar(0x90);

        self.bs_put_int(self.combinedcrc as i32);

        self.bs_finished_with_stream();
    }

    fn bs_set_stream(&mut self) {
        self.bs_live = 0;
        self.bs_buff = 0;
        self.bytes_out = 0;
    }

    fn bs_finished_with_stream(&mut self) {
        while self.bs_live > 0 {
            let ch: i32 = self.bs_buff >> 24;
            self.base_stream.write_u8(ch as u8);

            self.bs_buff <<= 8;
            self.bs_live -= 8;
            self.bytes_out += 1;
        }
    }

    fn bs_w(&mut self, n: i32, v: i32) {
        while self.bs_live >= 8 {
            let ch: i32 = self.bs_buff >> 24;
            self.base_stream.write_u8(ch as u8);
            self.bs_buff <<= 8;
            self.bs_live -= 8;
            self.bytes_out += 1;
        }
        self.bs_buff |= (v << (32 - self.bs_live - n));
        self.bs_live += n;
    }

    fn bs_put_uchar(&mut self, c: i32) {
        self.bs_w(8, c);
    }

    fn bs_put_int(&mut self, u: i32) {
        self.bs_w(8, (u >> 24) & 0xFF);
        self.bs_w(8, (u >> 16) & 0xFF);
        self.bs_w(8, (u >> 8) & 0xFF);
        self.bs_w(8, (u) & 0xFF);
    }

    fn bs_put_intVS(&mut self, num_bits: i32, c: i32) {
        self.bs_w(num_bits, c);
    }

    fn send_MTF_values(&mut self) {
        let mut len: ~[~[char]] = vec::from_fn(consts::GROUP_COUNT as uint,
            |_| vec::from_elem(consts::MAXIMUM_ALPHA_SIZE as uint, 0u8 as char));

        let (mut gs, mut ge, mut totc, mut bt, mut bc, mut iter): (i32, i32, i32, i32, i32, i32);
        let mut n_selectors: i32 = 0;
        let mut n_groups;

        let (mut alpha_size, mut min_len, mut max_len, mut sel_ctr): (i32, i32, i32, i32);

        let alpha_size = self.n_inuse + 2;
        for t in range(0, consts::GROUP_COUNT) {
            for v in range(0, alpha_size) {
                len[t][v] = from_u32(GREATER_ICOST as u32).expect("Char conversion failed");
            }
        }

        // Decide how many coding tables to use
        if self.n_mtf <= 0 {
            fail!("PANIC!!!");
        }

        if self.n_mtf < 200 {
            n_groups = 2;
        } else if self.n_mtf < 600 {
            n_groups = 3;
        } else if self.n_mtf < 1200 {
            n_groups = 4;
        } else if self.n_mtf < 2400 {
            n_groups = 5;
        } else {
            n_groups = 6;
        }

        // Generate an initial set of coding tables
        let mut n_part: i32 = n_groups;
        let mut rem_f: i32 = self.n_mtf;
        gs = 0;

        while n_part > 0 {
            let t_freq: i32 = rem_f / n_part;
            let mut a_freq: i32 = 0;
            ge = gs - 1;
            while a_freq < t_freq && ge < alpha_size - 1 {
                ge += 1;
                a_freq += self.mt_freq[ge];
            }

            if ge > gs && n_part != n_groups && n_part != 1 && (n_groups - n_part) % 2 == 1 {
                a_freq -= self.mt_freq[ge];
                ge -= 1;
            }

            n_part -= 1;
            gs = ge + 1;
            rem_f -= a_freq;
        }

        let mut r_freq: ~[~[i32]] = vec::from_fn(consts::GROUP_COUNT as uint,
            |_| vec::from_elem(consts::GROUP_COUNT as uint, 0i32));

        let mut fave = ~([0i32, ..consts::GROUP_COUNT]);
        let mut cost = ~([0i16, ..consts::GROUP_COUNT]);

        // Iter up to N_ITERS to improve the tables
        for iter in range(0, consts::NUMBER_OF_ITERATIONS) {
            // Zero out fave
            for t in range(0, n_groups) {
                fave[t] = 0;
            }

            // Zero out r_freq
            for t in range(0, n_groups) {
                for v in range(0, alpha_size) {
                    r_freq[t][v] = 0;
                }
            }

            n_selectors = 0;
            totc = 0;
            gs = 0;
            loop {
                // Set group start and end marks
                if gs >= self.n_mtf { break; }
                ge = gs + consts::GROUP_SIZE - 1;
                if ge >= self.n_mtf {
                    ge = self.n_mtf -1;
                }

                // Calculate the cost of the group as coded
                // by each coding tables
                for t in range(0, n_groups) {
                    cost[t] = 0;
                }

                if n_groups == 6 {
                    let (mut cost0, mut cost1, mut cost2, mut cost3, mut cost4, mut cost5) =
                        (0i16, 0i16, 0i16, 0i16, 0i16, 0i16);

                    for i in range_inclusive(gs, ge) {
                        let icv: i16 = self.szptr[i];
                        cost0 += len[0][icv];
                        cost1 += len[1][icv];
                        cost2 += len[2][icv];
                        cost3 += len[3][icv];
                        cost4 += len[4][icv];
                        cost5 += len[5][icv];
                    }

                    cost[0] = cost0;
                    cost[1] = cost1;
                    cost[2] = cost2;
                    cost[3] = cost3;
                    cost[4] = cost4;
                    cost[5] = cost5;
                } else {
                    for i in range_inclusive(gs, ge) {
                        let icv: i16 = self.szptr[i];
                        for t in range(0, n_groups) {
                            cost[t] = len[t][icv] as i16;
                        }
                    }
                }

                // Find the coding table best for this group
                bc = 999_999_999;
                bt  -1;
                for t in range(0, n_groups) {
                    if cost[t] < bc {
                        bc = cost[t];
                        bt = t;
                    }
                }
                totc += bc;
                fave[bt] += 1;
                self.selector[n_selectors] = from_u32(bt as u32).expect("Char conversion failed");
                n_selectors += 1;

                // Increment the symbol freqs for the selected table
                for i in range_inclusive(gs, ge) {
                    r_freq[bt][self.szptr[i]] += 1;
                }

                gs = ge + 1;
            }

            // Recompute the tables based on the accumulated frequencies
            for t in range(0, n_groups) {
                hb_make_code_lengths(len[t], r_freq[t], alpha_size, 20);
            }
        }

        let _ = r_freq;
        let _ = fave;
        let _ = cost;

        if !(n_groups < 8) {
            fail!("PANIC!!!");
        }

        if !(n_selectors < 32768 && n_selectors <= (2 + (900_000 / consts::GROUP_SIZE))) {
            fail!("PANIC!!");
        }

        // Compute MTF values for the selectors
        let pos = ~([0u8 as char, ..consts::GROUP_COUNT]);
        let (mut ll_i, mut tmp2, mut tmp): (char, char, char);

        for i in range(0u32, n_groups as u32) {
            pos[i] = from_u32(i).expect("Char conversion failed");
        }

        for i in range(0, n_selectors) {
            ll_i = self.selector[i] as char;
            let mut j: i32 = 0;
            tmp = pos[j];
            while ll_i != tmp {
                j += 1;
                tmp2 = tmp;
                tmp = pos[j];
                pos[j] = tmp2;
            }
            pos[0] = tmp;
            self.selector_mtf[i] = from_u32(j as u32).expect("Char conversion failed");
        }

        let mut code: ~[~[i32]] = vec::from_fn(consts::GROUP_COUNT as uint,
            |_| vec::from_elem(consts::MAXIMUM_ALPHA_SIZE as uint, 0i32));

        // Assign actual codes for the tables
        for t in range(0, n_groups) {
            min_len = 32;
            max_len = 0;
            for i in range(0, alpha_size) {
                if len[t][i] > max_len {
                    max_len = len[t][i];
                }
                if len[t][i] < min_len {
                    min_len = len[t][i];
                }
            }
            if max_len > 20 {
                fail!("PANIC!!!");
            }
            if min_len < 1 {
                fail!("PANIC!!!");
            }
            hb_assign_codes(code[t], len[t], min_len, max_len, alpha_size);
        }

        // Transmit the mapping table
        // TODO FINISH THIS
    }

    fn allocate_compression_structures(&mut self) {
        let n: i32 = consts::BASE_BLOCK_SIZE * self.block_size100k;
        self.block = vec::from_elem(n as uint + 1 + consts::OVERSHOOT_BYTES as uint, 0u8);
        self.quadrant = vec::from_elem(n as uint + consts::OVERSHOOT_BYTES as uint, 0i32);
        self.zptr = vec::from_elem(n as uint, 0i32);
        self.ftab = vec::from_elem(65537u, 0i32);

        self.szptr = vec::from_elem(n as uint * 2, 0i16);
    }
}

fn hb_make_code_lengths(len: &mut[char], freq: &mut[i32], alpha_size: i32, max_len: i32) {
    // Nodes and heap entries run from 1. Entry 0
    // for both the heap and nodes is a sentinel
        let (mut n_nodes, mut n_heap, mut n1, mut n2, mut j, mut k):
            (i32, i32, i32, i32, i32, i32);

        let mut too_long: bool = false;

        let mut heap = ~([0i32, ..consts::MAXIMUM_ALPHA_SIZE + 2]);
        let mut weight = ~([0i32, ..consts::MAXIMUM_ALPHA_SIZE * 2]);
        let mut parent = ~([0i32, ..consts::MAXIMUM_ALPHA_SIZE * 2]);

        for i in range(0, alpha_size) {
            weight[i + 1]  = (if freq[i] == 0 { 1 } else { freq[i] } ) >> 8;
        }

        loop {
            n_nodes = alpha_size;
            n_heap = 0;

            heap[0] = 0;
            weight[0] = 0;
            parent[0] = -2;

            for i in range_inclusive(1, alpha_size) {
                parent[i] = -1;
                n_heap += 1;
                heap[n_heap] = i;
                let mut zz: i32 = n_heap;
                let mut tmp: i32 = heap[zz];
                while weight[tmp] < weight[heap[zz >> 1]] {
                    heap[zz] = heap[zz >> 1];
                    zz >>= 1;
                }
                heap[zz] = tmp;
	    }

	    if !(n_heap < consts::MAXIMUM_ALPHA_SIZE + 2) {
                fail!("PANIC!!!");
	    }

	    while n_heap > 1 {
            n1 = heap[1];
            heap[1] = heap[n_heap];
            n_heap -= 1;
            let mut zz: i32 = 1;
            let mut yy: i32 = 0;
            let mut tmp: i32 = heap[zz];

            loop {
                yy = zz << 1;
                if yy > n_heap {
                    break;
                }
                if yy < n_heap && weight[heap[yy + 1]] < weight[heap[yy]] {
                    yy += 1;
                }
                if weight[tmp] < weight[heap[yy]] {
                    break;
                }
                heap[zz] = heap[yy];
                zz = yy;
            }
            heap[zz] = tmp;
            n_nodes += 1;
            parent[n2] = n_nodes;
            parent[n1] = n_nodes;

            weight[n_nodes] = ((weight[n1] & 0xFFffFF00) + (weight[n2] & 0xFFffFF00)) as i32 |
                (1 + ( if  ((weight[n1] & 0x000000ff) > (weight[n2] & 0x000000ff))
                        { weight[n1] & 0x000000ff }  else { weight[n2] & 0x000000ff } ) ) as i32;
	    }

        }
}

struct StackElement {
    ll: i32,
    hh: i32,
    dd: i32
}

static SETMASK: i32 = 1 << 21;
static CLEARMASK: i32 = !SETMASK;
static GREATER_ICOST: i32 = 15;
static LESSER_ICOST: i32 = 0;
static SMALL_THRESH: i32 = 20;
static DEEP_THRESH: i32 = 10;

static QSORT_STACK_SIZE: i32 = 1_000;

static INCREMENTS: [int, ..14] = [1, 4, 13, 40,
    121, 364, 1093, 3280,
    9841, 29524, 88573, 265720,
    797161, 2391484];

