const FOO: uint = 32;

fn main() {
    let arr = [0i, ..FOO];
    for n in arr.iter() {
        println!("n: {}", n);
    }
}