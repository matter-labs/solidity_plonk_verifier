# Code Generator for Solidity Verifier
## Generate
` cargo test -- --nocapture test_render_verification_key` 

## Test
Enter hardhat director 

`cd hardhat/`

run test with default proof file

`hardhat test`

run test with geth dev mode

`hardhat test --network geth`

run test with custom proof file

`PROOF_FILE=./../custom.proof hardhat test --network geth`