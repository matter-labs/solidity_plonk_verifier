# Code Generator for Solidity Verifier
## Install code generator tool
```
git clone ssh://git@github.com/matter-labs/solidity_plonk_verifier
cd solidity_plonk_verifier
cargo build
```

## Options
```bash
./target/debug/solidity_plonk_verifier --help
solidity_plonk_verifier 0.1.0

USAGE:
    solidity_plonk_verifier [OPTIONS] --verification-key <verification-key>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --main-gate <main-gate>
            Type of main gate [std_width4 | selector_optimized_width4] [default: std_width4]

        --output <output>                        Output directory [default: ./hardhat/contracts]
        --verification-key <verification-key>    Path to verification key(required)
```

## Generate for Standard Width4 Main Gate

### Generate
`./target/debug/solidity_plonk_verifier --verification-key ./block_vk_20_keccak.key`

### Test
Enter hardhat directory

`cd hardhat/`

run test 

`PROOF_FILE=./block_proof_20_keccak.proof hardhat test`

## Generate for Selector Optimized Width4 Main Gate

### Generate
`./target/debug/solidity_plonk_verifier --verification-key ./selector_optimized_vk_keccak.key --main-gate selector_optimized_width4`

### Test
Enter hardhat directory

`cd hardhat/`

run test 

`PROOF_FILE=./selector_optimized_proof_keccak.proof hardhat test`
