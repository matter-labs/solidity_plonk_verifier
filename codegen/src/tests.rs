use franklin_crypto::{bellman::{Engine, Field, PrimeField, SynthesisError, bn256::{Bn256, Fr}, groth16::create_proof, kate_commitment::{Crs, CrsForMonomialForm}, plonk::{better_better_cs::{                
                cs::{
                    Assembly, Circuit, ConstraintSystem, Gate, GateInternal,
                    LookupTableApplication, PlonkCsWidth4WithNextStepAndCustomGatesParams,
                    PolyIdentifier, Setup, TrivialAssembly,
                },
                setup::VerificationKey,
                gates::selector_optimized_with_d_next::SelectorOptimizedWidth4MainGateWithDNext,
                proof::Proof,
                verifier,
            }, commitments::transcript::{Transcript, keccak_transcript::RollingKeccakTranscript}}, worker::Worker}, plonk::circuit::{
        allocated_num::{AllocatedNum, Num},
        boolean::Boolean,
        custom_rescue_gate::Rescue5CustomGate,
        tables::inscribe_default_range_table_for_bit_width_over_first_three_columns,
    }};

use crate::{circuits::{SelectorOptimizedDummyCircuit, DummyCircuit}, generate};

fn generate_setup_vk_and_proof<E: Engine, C: Circuit<E>, T: Transcript<E::Fr>>(
    circuit: &C, 
    transcript_params: Option<T::InitializationParameters>,
) -> Result<(), SynthesisError> {
    let mut cs = TrivialAssembly::<
        E,
        PlonkCsWidth4WithNextStepAndCustomGatesParams,
        SelectorOptimizedWidth4MainGateWithDNext,
    >::new();

    circuit.synthesize(&mut cs)?;
    cs.finalize();

    assert!(cs.is_satisfied());

    let worker = Worker::new();

    let setup = cs.create_setup::<C>(&worker)?;

    println!("domain size {}", setup.n);

    let domain_size = setup.n.clone();

    let crs = Crs::<E, CrsForMonomialForm>::crs_42(domain_size.next_power_of_two(), &worker);

    let vk = VerificationKey::from_setup(&setup, &worker, &crs)?;
    println!("{:#?}", vk);

    let vk_file_name = "/tmp/selector_optimized_vk_keccak.key";
    let proof_file_name = "/tmp/selector_optimized_proof_keccak.proof";

    let mut vk_writer = std::fs::File::create(vk_file_name).expect("create vk file");
    vk.write(&mut vk_writer).expect("write vk into file");

    let proof = cs.create_proof::<_, T>(&worker, &setup, &crs, transcript_params.clone())?;

    let mut proof_writer = std::fs::File::create(proof_file_name).expect("create proof file");
    proof
        .write(&mut proof_writer)
        .expect("write proof into file");

    let verified = verifier::verify::<E, _, T>(&vk, &proof, transcript_params)?;

    assert!(verified, "proof verification failed");

    Ok(())
}

#[test]
fn test_create_vk_and_proof_file_for_std_main_gate() {
    let circuit = DummyCircuit;
    generate_setup_vk_and_proof::<Bn256, _, RollingKeccakTranscript<<Bn256 as Engine>::Fr>>(&circuit).unwrap()
}

#[test]
fn test_create_vk_and_proof_file_for_selector_optimized_main_gate() {
    let circuit = SelectorOptimizedDummyCircuit;

    generate_setup_vk_and_proof::<Bn256, _, RollingKeccakTranscript<<Bn256 as Engine>::Fr>>(&circuit).unwrap()
}

#[test]
fn test_render_vk_with_default_main_gate() {
    generate(
        "./../block_vk_20_keccak.key",
        "./../hardhat/contracts",
        None,
    );
}

#[test]
fn test_generate_vk_and_proof_for_dummy_circuit_with_selector_optimized_main_gate() {
    generate_for_selector_optimized(
        "./../block_vk_20_keccak.key",
        "./../hardhat/contracts",
        None,
    );
}
