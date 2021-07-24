pragma solidity >0.7;

import "./Pairing.sol";

uint256 constant STATE_WIDTH = 4;
uint256 constant NUM_SELECTORS = 7;
uint256 constant NUM_GATES = 2;

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
        vk.domain_size = 512;
        vk.omega = PairingsBn254.new_fr(0xfc88fce4f47faf73a154901e28615b16ac07c89b02332a5d5573c1d89a0bd30d);
        vk.gate_setup_commitments[0] = PairingsBn254.new_g1(0x4e11560bccfab22f1046098a83eb48938cbe51c75da9a6a284bb296e419a9411,0x615f6fd5d46bfde46558ef15a9d4ef62c04a10ef6ecce13f1e459fce53221b1e);
        vk.gate_setup_commitments[1] = PairingsBn254.new_g1(0x4ff8ca56bfe90eea0770fbac11a97ec75dcc717b887de8bcbda886eca0a5d125,0x0002b7efb5e0513096bfa4dfbd87bd55af7929e8f17595fdf7d9b2ac217e090d);
        vk.gate_setup_commitments[2] = PairingsBn254.new_g1(0xde035bcf3253dce869db46f6d401590c00b055e56ae5b8f4d034cdba2add6b2a,0xae44b48718eb25e370c03c354d40c2b9ea6342e19fd0e89ee4fb4900da70581b);
        vk.gate_setup_commitments[3] = PairingsBn254.new_g1(0x28e73823d9d30251cb4622ce12bf3202a382b165c09955c657a1fa517869d728,0x15636019b77f92a99c2a4138199ccabb31568134fbae05e04e72c01dd8c03501);
        vk.gate_setup_commitments[4] = PairingsBn254.new_g1(0xc5771a86a5c4cb489b933eefdc23fab3a2aeb36f7721bda540b9c465e25e1b2f,0x04b23ed6160cafc2d03f84e353c655f0207da3bb65fd3fd97580f35c15671329);
        vk.gate_setup_commitments[5] = PairingsBn254.new_g1(0x7a2ae6c9eee441b8b46c9e8d45a714ee5b303ddbb6ae0a70a14ae50573cdef0f,0x7bed4676a038ea159f35aa54a35620eb8673d55a409e0324e1796e62c3fe7923);
        vk.gate_setup_commitments[6] = PairingsBn254.new_g1(0x4d782194b76f52692e89afc33b5005fb4220c394d828fb68c31bbc4b5008b523,0xee439000a0ed6a91ce2c0769dbc307dd0308973b67c5462dfe90130c0e978702);
    }

}