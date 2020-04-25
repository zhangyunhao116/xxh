use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::SystemTime;

const CAP: usize = 256 * 1024;

fn main() {
    let start_time = SystemTime::now();
    let args: Vec<String> = env::args().collect();

    let filename = &args[1];
    let file = File::open(filename).expect("Invalid file path");

    let mut reader = BufReader::with_capacity(CAP, file);
    let mut digest = xxh::Xxh64::default();
    loop {
        let length = {
            let data = reader.fill_buf();
            match data {
                Err(_) => break,
                Ok(data) => {
                    digest.write(data);
                    data.len()
                },
            }
        };
        if length == 0 {
            break
        }
        reader.consume(length)
    }
    let result = digest.finish();
    println!("Finished `{}` in {}s\r\n\
    DEC: {}\r\n\
    HEX: {:x}", filename, SystemTime::now().duration_since(start_time).expect("Invalid system time").as_secs_f32(),
             result, result);
}
