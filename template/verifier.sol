pragma solidity >0.7;

import "./Pairing.sol";
import "./Verifier.sol";

uint256 constant STATE_WIDTH = {{STATE_WIDTH}};
uint256 constant NUM_SELECTORS = {{NUM_SELECTORS}};
uint256 constant NUM_GATES = {{NUM_GATES}};
uint256 constant NUM_G2_ELS = {{NUM_G2_ELS}};
uint256 constant NUM_LOOKUP_TABLES = {{NUM_LOOKUP_TABLES}};
uint256 constant SERIALIZED_PROOF_LENGTH = {{SERIALIZED_PROOF_LENGTH}};

struct VerificationKey {
    uint256 domain_size;
    uint256 num_inputs;
    PairingsBn254.Fr omega;
    PairingsBn254.G1Point[NUM_GATES] gate_selectors_commitments;
    PairingsBn254.G1Point[NUM_SELECTORS] gate_setup_commitments; // STATE_WIDTH for witness + multiplication + constant
    PairingsBn254.G1Point[STATE_WIDTH] permutation_commitments;
    uint256 total_lookup_entries_length;
    PairingsBn254.G1Point lookup_selector_commitment;
    PairingsBn254.G1Point[NUM_LOOKUP_TABLES] lookup_tables_commitments;
    PairingsBn254.G1Point lookup_table_type_commitment;

    PairingsBn254.Fr[STATE_WIDTH-1] non_residues;
    PairingsBn254.G2Point[NUM_G2_ELS] g2_elements;
}

contract Verifier is Plonk4VerifierWithAccessToDNext{

    function get_verification_key() internal pure returns(VerificationKey memory vk){
        vk.domain_size = {{domain_size}};
        vk.omega = {{domain_generator.el}};
        // coefficients
        {{#each gate_setup_commitments}}
        vk.gate_setup_commitments[{{@index}}] = {{this.g1}};
        {{/each}}
        // gate selectors
        {{#each gate_selectors_commitments}}
        vk.gate_selectors_commitments[{{@index}}] = {{this.g1}};
        {{/each}}
        // permutation
        {{#each permutation_commitments}}
        vk.permutation_commitments[{{@index}}] = {{this.g1}};
        {{/each}}
        vk.total_lookup_entries_length = {{total_lookup_entries_length}};
        
        // lookup table commitments
        vk.total_lookup_entries_length = {{total_lookup_entries_length}};
        vk.lookup_selector_commitment = {{lookup_selector_commitment.g1}};
        {{#each lookup_tables_commitments}}
        vk.lookup_tables_commitments[{{@index}}] = {{this.g1}};
        {{/each}}
        vk.lookup_table_type_commitment = {{lookup_table_type_commitment.g1}};
        
        // non residues
        {{#each non_residues}}
        vk.non_residues[{{@index}}] = {{this.el}};
        {{/each}}
        
        // g2 elements
        {{#each g2_elements}}
        vk.g2_elements[{{@index}}] = {{this.g2}};
        {{/each}}
    }


    function deserialize_proof(
        uint256[] memory public_inputs, 
        uint256[] memory serialized_proof
    ) internal pure returns(Proof memory proof) {
        require(serialized_proof.length == SERIALIZED_PROOF_LENGTH);
        proof.input_values = new uint256[](public_inputs.length);
        for (uint256 i = 0; i < public_inputs.length; i++) {
            proof.input_values[i] = public_inputs[i];
        }
 
        uint256 j = 0;
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.state_polys_commitments[i] = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
            );

            j += 2;
        }
        
        proof.copy_permutation_grand_product_commitment = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;
        
        proof.lookup_s_poly_commitment = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;
        
        proof.lookup_grand_product_commitment = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;
        
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.quotient_poly_parts_commitments[i] = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
            );

            j += 2;
        }
        
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.state_polys_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }
        
        for (uint256 i = 0; i < proof.state_polys_openings_at_z_omega.length; i++) {
            proof.state_polys_openings_at_z_omega[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }

        for (uint256 i = 0; i < proof.gate_setup_openings_at_z.length; i++) {
            proof.gate_setup_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }
        for (uint256 i = 0; i < proof.gate_selectors_openings_at_z.length; i++) {
            proof.gate_selectors_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }
        for (uint256 i = 0; i < proof.copy_permutation_polys_openings_at_z.length; i++) {
            proof.copy_permutation_polys_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }
        
        proof.copy_permutation_grand_product_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j += 1;

        proof.lookup_s_poly_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j += 1;
        proof.lookup_grand_product_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j += 1;

        proof.lookup_t_poly_opening_at_z = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j += 1;

        proof.lookup_t_poly_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j += 1;
        proof.lookup_selector_poly_opening_at_z = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j += 1;
        proof.lookup_table_type_poly_opening_at_z = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j += 1;

        proof.quotient_poly_opening_at_z = PairingsBn254.new_fr(
            serialized_proof[j]
        );

        j += 1;

        proof.linearization_poly_opening_at_z = PairingsBn254.new_fr(
            serialized_proof[j]
        );

        j += 1;

        proof.opening_proof_at_z = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;

        proof.opening_proof_at_z_omega = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
        );
    }
    
    function verify_serialized_proof(
        uint256[] memory public_inputs, 
        uint256[] memory serialized_proof
    ) public view returns (bool) {
        VerificationKey memory vk = get_verification_key();
        require(vk.num_inputs == public_inputs.length);

        Proof memory proof = deserialize_proof(public_inputs, serialized_proof);

        bool valid = verify(proof, vk);

        return valid;
    }  
}