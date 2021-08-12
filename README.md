# Code Generator for Solidity Verifier
## Install code generator tool
```
git clone ssh://git@github.com/matter-labs/solidity_plonk_verifier
cd solidity_plonk_verifier
```

## Options
```bash
cargo run -p codegen-bin generate -- --help
solidity_plonk_verifier 0.1.0

USAGE:
    generate [OPTIONS] --verification-key <verification-key>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --output <output>                        Output directory [default: ./hardhat/contracts]
        --verification-key <verification-key>    Path to verification key(required)
```

## Generate
`cargo run -p codegen-bin generate  -- --verification-key ./data/block_vk_20_keccak.key`

### Test
Enter hardhat directory

`cd hardhat/`

run test 

`PROOF_FILE=./block_proof_20_keccak.proof hardhat test`
