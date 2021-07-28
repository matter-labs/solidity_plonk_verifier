pragma solidity >0.7;

import "./Pairing.sol";
import "./Verifier.sol";

uint256 constant STATE_WIDTH = 4;
uint256 constant NUM_SELECTORS = 7;
uint256 constant NUM_GATES = 2;
uint256 constant NUM_G2_ELS = 2;
uint256 constant NUM_LOOKUP_TABLES = 4;
uint256 constant SERIALIZED_PROOF_LENGTH = 0;

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
        vk.domain_size = 512;
        vk.omega = PairingsBn254.new_fr(0xfc88fce4f47faf73a154901e28615b16ac07c89b02332a5d5573c1d89a0bd30d);
        // coefficients
        vk.gate_setup_commitments[0] = PairingsBn254.new_g1(0x4e11560bccfab22f1046098a83eb48938cbe51c75da9a6a284bb296e419a9411,0x615f6fd5d46bfde46558ef15a9d4ef62c04a10ef6ecce13f1e459fce53221b1e);
        vk.gate_setup_commitments[1] = PairingsBn254.new_g1(0x4ff8ca56bfe90eea0770fbac11a97ec75dcc717b887de8bcbda886eca0a5d125,0x0002b7efb5e0513096bfa4dfbd87bd55af7929e8f17595fdf7d9b2ac217e090d);
        vk.gate_setup_commitments[2] = PairingsBn254.new_g1(0xde035bcf3253dce869db46f6d401590c00b055e56ae5b8f4d034cdba2add6b2a,0xae44b48718eb25e370c03c354d40c2b9ea6342e19fd0e89ee4fb4900da70581b);
        vk.gate_setup_commitments[3] = PairingsBn254.new_g1(0x28e73823d9d30251cb4622ce12bf3202a382b165c09955c657a1fa517869d728,0x15636019b77f92a99c2a4138199ccabb31568134fbae05e04e72c01dd8c03501);
        vk.gate_setup_commitments[4] = PairingsBn254.new_g1(0xc5771a86a5c4cb489b933eefdc23fab3a2aeb36f7721bda540b9c465e25e1b2f,0x04b23ed6160cafc2d03f84e353c655f0207da3bb65fd3fd97580f35c15671329);
        vk.gate_setup_commitments[5] = PairingsBn254.new_g1(0x7a2ae6c9eee441b8b46c9e8d45a714ee5b303ddbb6ae0a70a14ae50573cdef0f,0x7bed4676a038ea159f35aa54a35620eb8673d55a409e0324e1796e62c3fe7923);
        vk.gate_setup_commitments[6] = PairingsBn254.new_g1(0x4d782194b76f52692e89afc33b5005fb4220c394d828fb68c31bbc4b5008b523,0xee439000a0ed6a91ce2c0769dbc307dd0308973b67c5462dfe90130c0e978702);
        // gate selectors
        vk.gate_selectors_commitments[0] = PairingsBn254.new_g1(0xa41b032d545d11e092401220ce9aab9f0cd4c4830b28cc7b81f06c44c14fbb0f,0x2eccf9d0e0f96882f9458b7438d2e848ce5dfce6023a7d86a82d6ef0bedeba03);
        vk.gate_selectors_commitments[1] = PairingsBn254.new_g1(0x883b4b4ad120ca659321dc25415a4f4d12ed19027e683e3465bc29cbbb0a9d11,0x68dc98e785ef0358b016f35735d549cec3165e398d47c4a45e3162823f8d8a23);
        // permutation
        vk.permutation_commitments[0] = PairingsBn254.new_g1(0x31bbc7639c756737d08e87b4c155e41621b7bd2d9ab0887c31eddfada468261d,0x6881aabdc75ddc7431116a23714501e84db99a4143fa1b74542df2d3e003eb28);
        vk.permutation_commitments[1] = PairingsBn254.new_g1(0x92b3d9262b5441d721896896dd81c902ed29aa67115d38381ca3032eee99be15,0xf943eefe5bb70fed6dfc2b41edb781070206fe50bdd5584308745074185d631e);
        vk.permutation_commitments[2] = PairingsBn254.new_g1(0x87196061c8074798a0902d9c8a045989b7284bada57736d7b906ad633aa02914,0x54dc989a2989abda9d38e6674ff992339191aca55775ce60d5bb1577b8afa415);
        vk.permutation_commitments[3] = PairingsBn254.new_g1(0x7e5bc2733b0052ba9b063dca5cb654f6997a5b50978e5eea34effadf770f7a0a,0x0a675512889a1de2b2fdd455479b855fc5dc376e3a7b0380c55ce3d20117bf28);
        vk.total_lookup_entries_length = 1;
        
        // lookup table commitments
        vk.total_lookup_entries_length = 1;
        vk.lookup_selector_commitment = PairingsBn254.new_g1(0xdf35ba05268d8cad71517ec7671d8337779f766f62378226236714d188df5f24,0xb34d011713f3f0e28c8e382fddba45d7ce1741ac8f4ae9d8742097f9e0a53502);
        vk.lookup_tables_commitments[0] = PairingsBn254.new_g1(0x02a25fa8bdf6a73a7bfc8c7286db5244bc7962f7818c66ce247d6fc6e99a771c,0x739203908f3191ad08e31576ef046177db7cf87870c6abf7c655adf935b16000);
        vk.lookup_tables_commitments[1] = PairingsBn254.new_g1(0x984304a57451362df6dd8ab8371ef12c5c94ba6acf8b94d28f15d2ecd5422b13,0xd6179d4dd47a01f60d619a970fa4718d6b5859852a5d1ed18472e6def20a2101);
        vk.lookup_tables_commitments[2] = PairingsBn254.new_g1(0xfafd089e793e656673ee139726d3c0b9ed50a95408ed03845925b6d6e51d890c,0x38fec2b504e6390efacbc4433296022a92939b36d598902b815ac149ae93c11a);
        vk.lookup_tables_commitments[3] = PairingsBn254.new_g1(0x13d33d7490b0e96bee7fd1da76fc26c4d99f796274ae6fb7cc2fbda08bd7b601,0x6a0d62001df418104f56ec2bc04f9d2b485af03e2552023d9ac50dd88cbdaf15);
        vk.lookup_table_type_commitment = PairingsBn254.new_g1(0x9b15e1ce6a3e053150c7f3e8633936debae55c5f89f43de6cf875db793ac1200,0x13032249af525cfa5d492c8e76d05bce1550710c7af9a0c4e26a648059758a1d);
        
        // non residues
        vk.non_residues[0] = PairingsBn254.new_fr(0x0500000000000000000000000000000000000000000000000000000000000000);
        vk.non_residues[1] = PairingsBn254.new_fr(0x0700000000000000000000000000000000000000000000000000000000000000);
        vk.non_residues[2] = PairingsBn254.new_fr(0x0a00000000000000000000000000000000000000000000000000000000000000);
        
        // g2 elements
        vk.g2_elements[0] = PairingsBn254.new_g2([0xedf692d95cbdde46ddda5ef7d422436779445c5e66006a42761e1f12efde0018,0xc212f3aeb785e49712e7a9353349aaf1255dfb31b7bf60723a480d9293938e19],[0xaa7dfa6601cce64c7bd3430c69e7d1e38f40cb8d8071ab4aeb6d8cdba55ec812,0x5b9722d1dcdaac55f38eb37033314bbc95330c69ad999eec75f05f58d0890609]);
        vk.g2_elements[1] = PairingsBn254.new_g2([0xe7e6af991b936dd662caae0472ba13801c5f3aa3ad44863d0f090d9ac8a86d11,0x919a5606733a9e0630a5e2d0b1670de93ce8cc6fb0496a7bb71596ba34097412],[0x72c1c48de52409d357618aac72fcbe14cf59f05162b5449630b6772e04416407,0xdb41a40b7cd0c6a84e832bd28bc179c9c5ac040dd0ed7d4a6bf8e516982d2225]);
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