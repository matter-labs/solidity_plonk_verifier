use codegen::generate;
use std::path::PathBuf;
use structopt::StructOpt;

const DEFAULT_OUTPUT_FILE: &str = "./hardhat/contracts";
const TEMPLATE_FILE_PATH: &str = "./codegen/template/verifier.sol";

#[derive(StructOpt, Debug)]
pub struct Opts {
    /// Path to verification key(required)
    #[structopt(long, parse(from_os_str))]
    verification_key: PathBuf,
    /// Output directory
    #[structopt(long, parse(from_os_str), default_value = DEFAULT_OUTPUT_FILE)]
    output: PathBuf,
}

fn main() {
    let opts = Opts::from_args();
    println!("{:#?}", opts);

    let Opts {
        verification_key,
        output,
    } = opts;

    generate(verification_key, output.clone(), Some(TEMPLATE_FILE_PATH));

    eprintln!("Success!");
}
