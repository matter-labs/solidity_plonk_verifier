pragma solidity >0.7;

import "./Pairing.sol";

uint256 constant STATE_WIDTH = {{STATE_WIDTH}};
uint256 constant NUM_SELECTORS = {{NUM_SELECTORS}};
uint256 constant NUM_GATES = {{NUM_GATES}};

struct VerificationKey {
    uint256 domain_size;
    uint256 num_inputs;
    PairingsBn254.Fr omega;
    PairingsBn254.G1Point[NUM_GATES] gate_selector_commitments;
    PairingsBn254.G1Point[NUM_SELECTORS] gate_setup_commitments; // STATE_WIDTH for witness + multiplication + constant
    PairingsBn254.G1Point[STATE_WIDTH] permutation_commitments;
    // PairingsBn254.G1Point[ACCESSIBLE_STATE_POLYS_ON_NEXT_STEP] next_step_gate_setup_commitments;
    uint256 total_lookup_entries_length;
    PairingsBn254.G1Point lookup_gate_setup_commitment;
    PairingsBn254.G1Point lookup_table_type_commitment;
    PairingsBn254.G1Point[] lookup_tables_commitments;

    PairingsBn254.Fr[STATE_WIDTH-1] permutation_non_residues;
    PairingsBn254.G2Point g2_x;
}

contract Verifier{

    function get_vk_key() internal pure returns(VerificationKey memory vk){
        vk.domain_size = {{domain_size}};
        vk.omega = {{domain_generator.fr}};
        {{#each setup_commitments}}
        vk.gate_setup_commitments[{{@index}}] = {{this.g1}};
        {{/each}}


    }

}