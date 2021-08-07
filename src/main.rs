use codegen::{generate, MainGateType};
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
    /// Type of main gate [std_width4 | selector_optimized_width4]
    #[structopt(long, default_value = "std_width4")]
    main_gate: MainGateType,
}

fn main() {
    let opts = Opts::from_args();    
    println!("{:#?}", opts);

    let Opts{
        verification_key,
        output,
        main_gate,
    } = opts;

    generate(main_gate, verification_key, output.clone(), Some(TEMPLATE_FILE_PATH));
    
    eprintln!("Success!");
}
