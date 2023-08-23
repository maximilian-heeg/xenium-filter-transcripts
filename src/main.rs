use clap::Parser;
use std::process;

use xenium_filter_transcripts::*;

fn main() {
    let args = Args::parse();

    if let Err(err) = run(args) {
        println!("{}", err);
        process::exit(1);
    }
}
