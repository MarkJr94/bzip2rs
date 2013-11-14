#[macro_escape];

macro_rules! do_while(
    ($body:expr, $cond:expr) => (
        loop {
            $body;

            if !($cond) { break; }
        }
    )
)