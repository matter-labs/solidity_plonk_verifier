// const { expect } = require("chai");
import * as chai from "chai";
import {ethers} from "hardhat";
import {Proof} from "./utils";
import {readFileSync, existsSync, PathLike} from "fs";
import { BigNumber, Contract } from "ethers";
import * as path from "path";

const expect = chai.expect;

const DEFAULT_PROOF_FILE = path.resolve(__dirname, "./../../block_proof_20_keccak.proof");

describe("Verifier", function () {

  let proof : Proof;
  let contract: Contract;

  before(async function(){
    let proofFile = process.env.PROOF_FILE || DEFAULT_PROOF_FILE;
    if(!existsSync(proofFile)){
      throw new Error(`proof file not found! ${proofFile}`);
      return;
    }    
    proof = Proof.read_from_file(proofFile);
    let factory = await ethers.getContractFactory("Verifier");
    contract = await factory.deploy();
    await contract.deployed();
  })

  it("Should verify proof", async function(){
    let [pub_inputs, encoding] = proof.serialize();
    let result = await contract.verify_serialized_proof(pub_inputs, encoding);
    expect(result, "proof verification failed").true;
  }); 
});
