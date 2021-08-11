#!/bin/bash
for main_gate in "std" "optimized"
do
    for vk_file in $(ls ./data/$main_gate/*.key)
    do
        proof_file=$(echo $vk_file | sed s/key/proof/g)
        proof_file=$(echo $proof_file | sed s/vk/proof/g)
        echo PROOF FILE $proof_file VK FILE $vk_file

        ./target/release/solidity_plonk_verifier --verification-key $vk_file
        cd hardhat
        PROOF_FILE=../$proof_file hh test
        cd -
    done
done