use codegen::generate;

const DEFAULT_OUTPUT_FILE: &str  = "./hardhat/contracts";
fn main() {
    let vk_file_path =std::env::args().nth(1).expect("should provide verification key path");
    let output_path = match std::env::args().nth(2){
        Some(path) => path.to_string(),
        None => DEFAULT_OUTPUT_FILE.to_string(),
    };

    let template_file_path = "./codegen/template/verifier.sol";

    generate(&vk_file_path, &output_path, Some(template_file_path)); 

    println!("Verifier saved into {}", output_path);
}
