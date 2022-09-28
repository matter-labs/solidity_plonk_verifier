# Code Generator for Solidity Verifier
## Install code generator tool
```
git clone ssh://git@github.com/matter-labs/solidity_plonk_verifier --branch dev
cd solidity_plonk_verifier
```

## Options
```bash
cargo run generate -- --help
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
`cargo run  -- --verification-key ./data/example.key`

### Test
Enter hardhat directory

`cd hardhat/`

run test 

`PROOF_FILE=../data/example.proof npm test`
or
`PROOF_FILE=../data/example.proof npx hardhat test`
