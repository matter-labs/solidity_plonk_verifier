pragma solidity ^0.8.0;

import "./Plonk4VerifierWithAccessToDNext.sol";
import "./UncheckedMath.sol";

contract Verifier is Plonk4VerifierWithAccessToDNext {
    using UncheckedMath for uint256;

    function get_verification_key() internal pure returns(VerificationKey memory vk){
        vk.num_inputs = 1;
        vk.domain_size = 512;
        vk.omega = PairingsBn254.new_fr(0x0dd30b9ad8c173555d2a33029bc807ac165b61281e9054a173af7ff4e4fc88fc);
        // coefficients
        vk.gate_setup_commitments[0] = PairingsBn254.new_g1(0x2c6005e59e769e19bb2f46bf2d9f65283a4b1c5c34aa2349862e4f1bb63ee91c,0x123ea7b388558f1d1b4c2213e9acff6476fafc7702f9a5b6855b312e34ddecec);
        vk.gate_setup_commitments[1] = PairingsBn254.new_g1(0x1bf67502fbb626a425575c04ba2d7e6001d0c79b78432f1618d830cdda31daf1,0x14b48b1ca4aa965a1b3b295f75db9dad001613638a2b8e94e1194ddc1dd3cdb6);
        vk.gate_setup_commitments[2] = PairingsBn254.new_g1(0x28a1c86f2bb7036a1e20477cb9bfb99e194344b229e0c1be5ca14d52d3d7d78b,0x2a870bd7c11012783e849c1ae6d24a2dac8ad38bb6ef4cbb25c31c001a2b6dcc);
        vk.gate_setup_commitments[3] = PairingsBn254.new_g1(0x2a81ec592f22c02e74a9885ff531b018274fafa7fa19dee1234f5191c5dce24b,0x27a3d4a6bfb3548d86de261c1a4d7a6826246cb9a5a17b252d70588cf9363510);
        vk.gate_setup_commitments[4] = PairingsBn254.new_g1(0x183d7488f5134ca26b0ce1e9fca43aeaf5457ecb9ca9c3bcd0751e21d0c10221,0x047c91abc16f5af01c6950aa55996c9929ad6bdf5b4228727495ee71eb2fc6c0);
        vk.gate_setup_commitments[5] = PairingsBn254.new_g1(0x1bba345553aca8b09d89bd9e534bb29ea6fc715a2f8c62fd5e3ffa6b64f7823b,0x2c26992dbb24abf84fae130447df736b372b93c46ac93f33a7b0c6d0477e1afa);
        vk.gate_setup_commitments[6] = PairingsBn254.new_g1(0x21b50c44947711561c80a88dd1bd2f12aa28f0d88df6d2060915179e74300fa3,0x0cd866091da9fd1bed85b8c85ba701ad01ed97142d1df1420312ae30090fe6c7);
        vk.gate_setup_commitments[7] = PairingsBn254.new_g1(0x0000000000000000000000000000000000000000000000000000000000000000,0x0000000000000000000000000000000000000000000000000000000000000001);
        // gate selectors
        vk.gate_selectors_commitments[0] = PairingsBn254.new_g1(0x1c4922a27df7b59b2b60afa704989db3b92793468239419f25058eeee54951b0,0x07566b3fb05764af536bb031d861eb5208f873e5622f9116afc209a70d521816);
        vk.gate_selectors_commitments[1] = PairingsBn254.new_g1(0x28ddac5a81f3045ba660e4c95962500814a130c77ac5f815d747e9a50ad125cf,0x2e6ad370b93c4de8b3211cf438f5a5808be895003d40d5346109b6d69d2a55c9);
        // permutation
        vk.permutation_commitments[0] = PairingsBn254.new_g1(0x08db50f68ace5a40979b76fd3c9e5fcfb2c769f33c2da8c068bfe142c51b74f0,0x03c4493634b2860db185269b3a67a0629d67a9f25e12e83702c649fa2b1c5be9);
        vk.permutation_commitments[1] = PairingsBn254.new_g1(0x1dbe246e6087c94e6525ffe99880530f55ae6c17c6157d07275f8dddbd9ab651,0x17d13a8038919d0f75eb7beb50dfaf1a3b1096954b7a3835a22b14484272dcf8);
        vk.permutation_commitments[2] = PairingsBn254.new_g1(0x218de33ba439143345ea15348362e6cbaefbd4e4de9ae3d99253f182ae0892e9,0x03a623176bfb207d053215c71078ba235a363c7607349ab9b8ba6060c30745b6);
        vk.permutation_commitments[3] = PairingsBn254.new_g1(0x1fb07394c8e800ed6f5193622a1574df742bb3f95e55d2e0dfcb5b80c33f3b3f,0x2e2ab2002bec17e39e49db8cafc5e3d0b16abce85ceb7332f69121b4bf7c30c6);
        // lookup table commitments
        vk.lookup_selector_commitment = PairingsBn254.new_g1(0x28a33ec6cc7ab335c0c3c742d845bb463132b5e9906b556e3a6500d61d7b4086,0x2de7bdeb5b21b1e11fdae54a178c71086944a0865555be0b7518f665b0466b89);
        vk.lookup_tables_commitments[0] = PairingsBn254.new_g1(0x1c779ae9c66f7d24ce668c81f76279bc4452db86728cfc7b3aa7f6bda85fa202,0x0060b135f9ad55c6f7abc67078f87cdb776104ef7615e308ad91318f90039273);
        vk.lookup_tables_commitments[1] = PairingsBn254.new_g1(0x132b42d5ecd2158fd2948bcf6aba945c2cf11e37b88addf62d365174a5044398,0x01210af2dee67284d11e5d2a8559586b8d71a40f979a610df6017ad44d9d17d6);
        vk.lookup_tables_commitments[2] = PairingsBn254.new_g1(0x0c891de5d6b625598403ed0854a950edb9c0d3269713ee7366653e799e08fdfa,0x1ac193ae49c15a812b9098d5369b93922a02963243c4cbfa0e39e604b5c2fe38);
        vk.lookup_tables_commitments[3] = PairingsBn254.new_g1(0x01b6d78ba0bd2fccb76fae7462799fd9c426fc76dad17fee6be9b090743dd313,0x15afbd8cd80dc59a3d0252253ef05a482b9d4fc02bec564f1018f41d00620d6a);
        vk.lookup_table_type_commitment = PairingsBn254.new_g1(0x1e1457a85814242ed4189182d7ca84cdcf8d07b091b2523dbdc19d2c0192c428,0x08e85e5a2c69c8616e18a0cf411837e54e9c9c1c0afd8a86fccaaaf6910efee8);
        // non residues
        vk.non_residues[0] = PairingsBn254.new_fr(0x0000000000000000000000000000000000000000000000000000000000000005);
        vk.non_residues[1] = PairingsBn254.new_fr(0x0000000000000000000000000000000000000000000000000000000000000007);
        vk.non_residues[2] = PairingsBn254.new_fr(0x000000000000000000000000000000000000000000000000000000000000000a);
        
        // g2 elements
        vk.g2_elements[0] = PairingsBn254.new_g2([0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2,0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed],[0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b,0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa]);
        vk.g2_elements[1] = PairingsBn254.new_g2([0x12740934ba9615b77b6a49b06fcce83ce90d67b1d0e2a530069e3a7306569a91,0x116da8c89a0d090f3d8644ada33a5f1c8013ba7204aeca62d66d931b99afe6e7],[0x25222d9816e5f86b4a7dedd00d04acc5c979c18bd22b834ea8c6d07c0ba441db,0x076441042e77b6309644b56251f059cf14befc72ac8a6157d30924e58dc4c172]);
    }

    function deserialize_proof(
        uint256[] calldata public_inputs, 
        uint256[] calldata serialized_proof
    ) internal pure returns(Proof memory proof) {
        // require(serialized_proof.length == 44); TODO
        proof.input_values = new uint256[](public_inputs.length);
        for (uint256 i = 0; i < public_inputs.length; i = i.uncheckedInc()) {
            proof.input_values[i] = public_inputs[i];
        }
 
        uint256 j;
        for (uint256 i = 0; i < STATE_WIDTH; i = i.uncheckedInc()) {
            proof.state_polys_commitments[i] = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
            );

            j = j.uncheckedAdd(2);
        }
        proof.copy_permutation_grand_product_commitment = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
        );
        j = j.uncheckedAdd(2);
        
        proof.lookup_s_poly_commitment = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
        );
        j = j.uncheckedAdd(2);

        proof.lookup_grand_product_commitment = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
        );
        j = j.uncheckedAdd(2);
        for (uint256 i = 0; i < proof.quotient_poly_parts_commitments.length; i = i.uncheckedInc()) {
            proof.quotient_poly_parts_commitments[i] = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
            );
            j = j.uncheckedAdd(2);
        }

        for (uint256 i = 0; i < proof.state_polys_openings_at_z.length; i = i.uncheckedInc()) {
            proof.state_polys_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j = j.uncheckedInc();
        }

        for (uint256 i = 0; i < proof.state_polys_openings_at_z_omega.length; i = i.uncheckedInc()) {
            proof.state_polys_openings_at_z_omega[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j = j.uncheckedInc();
        } 
        for (uint256 i = 0; i < proof.gate_selectors_openings_at_z.length; i = i.uncheckedInc()) {
            proof.gate_selectors_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j = j.uncheckedInc();
        }        
        for (uint256 i = 0; i < proof.copy_permutation_polys_openings_at_z.length; i = i.uncheckedInc()) {
            proof.copy_permutation_polys_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j = j.uncheckedInc();
        }
        proof.copy_permutation_grand_product_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j = j.uncheckedInc();
        proof.lookup_s_poly_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );
        j = j.uncheckedInc();
        proof.lookup_grand_product_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j = j.uncheckedInc();
        proof.lookup_t_poly_opening_at_z = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j = j.uncheckedInc();
        proof.lookup_t_poly_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );
        j = j.uncheckedInc();
        proof.lookup_selector_poly_opening_at_z = PairingsBn254.new_fr(
                serialized_proof[j]
            );
        j = j.uncheckedInc();
        proof.lookup_table_type_poly_opening_at_z = PairingsBn254.new_fr(
                serialized_proof[j]
            );
        j = j.uncheckedInc();
        proof.quotient_poly_opening_at_z = PairingsBn254.new_fr(
            serialized_proof[j]
        );
        j = j.uncheckedInc();
        proof.linearization_poly_opening_at_z = PairingsBn254.new_fr(
            serialized_proof[j]
        );
        j = j.uncheckedInc();
        proof.opening_proof_at_z = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
        );
        j = j.uncheckedAdd(2);
        proof.opening_proof_at_z_omega = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
        );
    }
    
    function verify_serialized_proof(
        uint256[] calldata public_inputs, 
        uint256[] calldata serialized_proof
    ) public view returns (bool) {
        VerificationKey memory vk = get_verification_key();
        require(vk.num_inputs == public_inputs.length);

        Proof memory proof = deserialize_proof(public_inputs, serialized_proof);

        return verify(proof, vk);
    }
}
