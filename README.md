# Code Generator for Solidity Verifier
## Generate

`cargo run -- ./block_vk_20_keccak.key`

Write solidity file into custom directory

`cargo run -- ./block_vk_20_keccak.key ./path/to/directory`

## Test
Enter hardhat director 

`cd hardhat/`

run test with default proof file

`hardhat test`

run test with geth dev mode

`hardhat test --network geth`

run test with custom proof file

`PROOF_FILE=./../custom.proof hardhat test --network geth`