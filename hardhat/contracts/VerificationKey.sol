pragma solidity >0.7;

uint256 constant STATE_WIDTH = 4;
uint256 constant NUM_G2_ELS = 2;

struct VerificationKey {
    uint256 domain_size;
    uint256 num_inputs;
    PairingsBn254.Fr omega;
    PairingsBn254.G1Point[2] gate_selectors_commitments;
    PairingsBn254.G1Point[8] gate_setup_commitments;
    PairingsBn254.G1Point[STATE_WIDTH] permutation_commitments;
    PairingsBn254.Fr[STATE_WIDTH-1] non_residues;
    PairingsBn254.G2Point[NUM_G2_ELS] g2_elements;
}

contract Plonk4VerifierWithAccessToDNext {
    using PairingsBn254 for PairingsBn254.G1Point;
    using PairingsBn254 for PairingsBn254.G2Point;
    using PairingsBn254 for PairingsBn254.Fr;
    
    using TranscriptLibrary for TranscriptLibrary.Transcript;
    
    struct Proof {
        uint256[] input_values;
        // commitments
        PairingsBn254.G1Point[STATE_WIDTH] state_polys_commitments;
        PairingsBn254.G1Point copy_permutation_grand_product_commitment;
        PairingsBn254.G1Point[STATE_WIDTH] quotient_poly_parts_commitments;
        
        // openings
        PairingsBn254.Fr[STATE_WIDTH] state_polys_openings_at_z;
        PairingsBn254.Fr[1] state_polys_openings_at_z_omega; // TODO: not use array while there is only D_next
        PairingsBn254.Fr[1] gate_selectors_openings_at_z;
        PairingsBn254.Fr[STATE_WIDTH-1] copy_permutation_polys_openings_at_z;
        PairingsBn254.Fr copy_permutation_grand_product_opening_at_z_omega;
        PairingsBn254.Fr quotient_poly_opening_at_z;
        PairingsBn254.Fr linearization_poly_opening_at_z;

        PairingsBn254.G1Point opening_proof_at_z;
        PairingsBn254.G1Point opening_proof_at_z_omega;
    }
    
    struct PartialVerifierState {
        PairingsBn254.Fr zero;
        PairingsBn254.Fr alpha;
        PairingsBn254.Fr beta;
        PairingsBn254.Fr gamma;
        PairingsBn254.Fr[6] alpha_values;
        PairingsBn254.Fr v;
        PairingsBn254.Fr u;
        PairingsBn254.Fr z;        
        PairingsBn254.Fr z_omega;
        PairingsBn254.Fr z_minus_last_omega;
        PairingsBn254.Fr l_0_at_z;
        PairingsBn254.Fr l_n_minus_one_at_z;
        PairingsBn254.Fr t;
        PairingsBn254.G1Point tp;
    }

    function evaluate_l0_at_point(
        uint256 domain_size, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr memory num) {
        PairingsBn254.Fr memory one = PairingsBn254.new_fr(1);

        PairingsBn254.Fr memory size_fe = PairingsBn254.new_fr(domain_size);
        PairingsBn254.Fr memory den = at.copy();
        den.sub_assign(one);
        den.mul_assign(size_fe);

        den = den.inverse();

        num = at.pow(domain_size);
        num.sub_assign(one);
        num.mul_assign(den);
    }
    
    function evaluate_lagrange_poly_out_of_domain(
        uint256 poly_num, 
        uint256 domain_size, 
        PairingsBn254.Fr memory omega, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr memory res) {
        // (omega^i / N) / (X - omega^i) * (X^N - 1)
        require(poly_num < domain_size);
        PairingsBn254.Fr memory one = PairingsBn254.new_fr(1);
        PairingsBn254.Fr memory omega_power = omega.pow(poly_num);
        res = at.pow(domain_size);
        res.sub_assign(one);
        require(res.value != 0); // Vanishing polynomial can not be zero at point `at`
        res.mul_assign(omega_power);
        
        PairingsBn254.Fr memory den = PairingsBn254.copy(at);
        den.sub_assign(omega_power);
        den.mul_assign(PairingsBn254.new_fr(domain_size));
        
        den = den.inverse();
        
        res.mul_assign(den);
    }    
    
    function evaluate_vanishing(
        uint256 domain_size, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr memory res) {
        res = at.pow(domain_size);
        res.sub_assign(PairingsBn254.new_fr(1));
    }
        
    function initialize_transcript(Proof memory proof, VerificationKey memory vk) internal pure returns (PartialVerifierState memory state) {
        TranscriptLibrary.Transcript memory transcript = TranscriptLibrary.new_transcript();

        for(uint256 i =0; i < vk.num_inputs; i++){
            transcript.update_with_u256(proof.input_values[i]);
        }

        for(uint256 i =0; i < STATE_WIDTH; i++){
            transcript.update_with_g1(proof.state_polys_commitments[i]);
        }

        state.beta = transcript.get_challenge();
        state.gamma = transcript.get_challenge();

        transcript.update_with_g1(proof.copy_permutation_grand_product_commitment);
        state.alpha = transcript.get_challenge();

        for(uint256 i =0; i < proof.quotient_poly_parts_commitments.length; i++){
            transcript.update_with_g1(proof.quotient_poly_parts_commitments[i]);
        }   
        state.z = transcript.get_challenge();

        transcript.update_with_fr(proof.quotient_poly_opening_at_z);

        for(uint256 i =0; i < proof.state_polys_openings_at_z.length; i++){
            transcript.update_with_fr(proof.state_polys_openings_at_z[i]);
        }        

        for(uint256 i =0; i < proof.state_polys_openings_at_z_omega.length; i++){
            transcript.update_with_fr(proof.state_polys_openings_at_z_omega[i]);
        }             
        for(uint256 i =0; i < proof.gate_selectors_openings_at_z.length; i++){
            transcript.update_with_fr(proof.gate_selectors_openings_at_z[i]);
        }          
        for(uint256 i =0; i < proof.copy_permutation_polys_openings_at_z.length; i++){
            transcript.update_with_fr(proof.copy_permutation_polys_openings_at_z[i]);
        }          

        state.z_omega = state.z.copy();
        state.z_omega.mul_assign(vk.omega);

        transcript.update_with_fr(proof.copy_permutation_grand_product_opening_at_z_omega);

        transcript.update_with_fr(proof.linearization_poly_opening_at_z);

        state.v = transcript.get_challenge();
        
        transcript.update_with_g1(proof.opening_proof_at_z);
        transcript.update_with_g1(proof.opening_proof_at_z_omega);

        state.u = transcript.get_challenge();
    }

    // compute some powers of challenge alpha([alpha^1, .. alpha^8]) 
    function compute_powers_of_alpha(PartialVerifierState memory state) public pure{
        require(state.alpha.value != 0);
        state.alpha_values[0] = PairingsBn254.new_fr(1);
        state.alpha_values[1] = state.alpha.copy();
        PairingsBn254.Fr memory current_alpha = state.alpha.copy();
        for(uint256 i=2; i < state.alpha_values.length; i++){
            current_alpha.mul_assign(state.alpha);
            state.alpha_values[i] = current_alpha.copy();
        }
    }

    function verify(Proof memory proof, VerificationKey memory vk) internal view returns (bool) {        
        // we initialize all challenges beforehand, we can draw each challenge in its own place
        PartialVerifierState memory state = initialize_transcript(proof, vk);
        if(verify_quotient_evaluation(vk, proof, state)== false){
                return false;
        }
        require(proof.state_polys_openings_at_z_omega.length == 1); // TODO

        
        PairingsBn254.G1Point memory quotient_result =  proof.quotient_poly_parts_commitments[0].copy_g1();
        {
            // block scope
            PairingsBn254.Fr memory z_in_domain_size = state.z.pow(vk.domain_size);
            PairingsBn254.Fr memory current_z = z_in_domain_size.copy();
            PairingsBn254.G1Point memory tp;
            // start from i =1 
            for(uint256 i = 1; i < proof.quotient_poly_parts_commitments.length; i++){
                tp = proof.quotient_poly_parts_commitments[i].copy_g1();
                tp.point_mul_assign(current_z);
                quotient_result.point_add_assign(tp);

                current_z.mul_assign(z_in_domain_size);
            }
        }
        
        Queries memory queries = prepare_queries(vk, proof, state);
        queries.commitments_at_z[0] = quotient_result;
        queries.values_at_z[0] = proof.quotient_poly_opening_at_z;
        queries.commitments_at_z[1] = aggregated_linearization_commitment(vk, proof, state);
        queries.values_at_z[1] = proof.linearization_poly_opening_at_z;

        require(queries.commitments_at_z.length == queries.values_at_z.length);

        PairingsBn254.G1Point memory  aggregated_commitment_at_z = queries.commitments_at_z[0];
        
        PairingsBn254.Fr memory  aggregated_opening_at_z = queries.values_at_z[0];
        PairingsBn254.Fr memory  aggregation_challenge = PairingsBn254.new_fr(1);
        PairingsBn254.G1Point memory scaled;
        for(uint256 i = 1; i < queries.commitments_at_z.length; i++){
            aggregation_challenge.mul_assign(state.v);
            scaled = queries.commitments_at_z[i].point_mul(aggregation_challenge);
            aggregated_commitment_at_z.point_add_assign(scaled);

            state.t = queries.values_at_z[i];
            state.t.mul_assign(aggregation_challenge);
            aggregated_opening_at_z.add_assign(state.t);
        }

        aggregation_challenge.mul_assign(state.v);

        PairingsBn254.G1Point memory  aggregated_commitment_at_z_omega = queries.commitments_at_z_omega[0].point_mul(aggregation_challenge);
        PairingsBn254.Fr memory  aggregated_opening_at_z_omega = queries.values_at_z_omega[0];        
        aggregated_opening_at_z_omega.mul_assign(aggregation_challenge);
        for(uint256 i = 1; i < queries.commitments_at_z_omega.length; i++){
            aggregation_challenge.mul_assign(state.v);

            scaled = queries.commitments_at_z_omega[i].point_mul(aggregation_challenge);
            aggregated_commitment_at_z_omega.point_add_assign(scaled);

            state.t = queries.values_at_z_omega[i];
            state.t.mul_assign(aggregation_challenge);
            aggregated_opening_at_z_omega.add_assign(state.t);
        }

        return final_pairing(
            vk.g2_elements,
            proof, 
            state, 
            aggregated_commitment_at_z,
            aggregated_commitment_at_z_omega,
            aggregated_opening_at_z,            
            aggregated_opening_at_z_omega
        );

    }

    function verify_quotient_evaluation(VerificationKey memory vk,  Proof memory proof, PartialVerifierState memory state) internal view returns(bool){        
        uint256[] memory lagrange_poly_numbers = new uint256[](vk.num_inputs);
        for (uint256 i = 0; i < lagrange_poly_numbers.length; i++) {
            lagrange_poly_numbers[i] = i;
        }
        // require(vk.num_inputs > 0); // TODO
        
        PairingsBn254.Fr memory inputs_term = PairingsBn254.new_fr(0);
        for(uint256 i =0; i < vk.num_inputs; i++){
            // TODO we may use batched lagrange compputation            
            state.t = evaluate_lagrange_poly_out_of_domain(i, vk.domain_size, vk.omega, state.z);            
            state.t.mul_assign(PairingsBn254.new_fr(proof.input_values[i]));
            inputs_term.add_assign(state.t);
        } 
        inputs_term.mul_assign(proof.gate_selectors_openings_at_z[0]);         
        PairingsBn254.Fr memory result = proof.linearization_poly_opening_at_z.copy();
        result.add_assign(inputs_term);
        
        // compute powers of alpha 
        compute_powers_of_alpha(state);
        PairingsBn254.Fr memory factor = state.alpha_values[4].copy();
        factor.mul_assign(proof.copy_permutation_grand_product_opening_at_z_omega);

        // - alpha_0 * (a + perm(z) * beta + gamma)*()*(d + gamma) * z(z*omega)
        require(proof.copy_permutation_polys_openings_at_z.length == STATE_WIDTH-1);
        PairingsBn254.Fr memory t;  // TMP;
        for(uint256 i = 0; i < proof.copy_permutation_polys_openings_at_z.length; i++){
            t = proof.copy_permutation_polys_openings_at_z[i].copy();
            t.mul_assign(state.beta);
            t.add_assign(proof.state_polys_openings_at_z[i]);
            t.add_assign(state.gamma);

            factor.mul_assign(t);
        }

        t = proof.state_polys_openings_at_z[3].copy();
        t.add_assign(state.gamma);
        factor.mul_assign(t);
        result.sub_assign(factor);
        
        // - L_0(z) * alpha_1
        PairingsBn254.Fr memory l_0_at_z = evaluate_l0_at_point(vk.domain_size, state.z);
        l_0_at_z.mul_assign(state.alpha_values[4 + 1]);
        result.sub_assign(l_0_at_z);

        PairingsBn254.Fr memory lhs = proof.quotient_poly_opening_at_z.copy();
        lhs.mul_assign(evaluate_vanishing(vk.domain_size, state.z));    
        return lhs.value == result.value;
    }
    function aggregated_linearization_commitment(VerificationKey memory vk,  Proof memory proof, PartialVerifierState memory state) internal view returns(PairingsBn254.G1Point memory result){                
        // qMain*(Q_a * A + Q_b * B + Q_c * C + Q_d * D + Q_m * A*B + Q_const + Q_dNext * D_next)
        result = PairingsBn254.new_g1(0, 0);
        // Q_a * A        
        PairingsBn254.G1Point memory scaled = vk.gate_setup_commitments[0].point_mul(proof.state_polys_openings_at_z[0]);
        result.point_add_assign(scaled);
        // Q_b * B
        scaled = vk.gate_setup_commitments[1].point_mul(proof.state_polys_openings_at_z[1]);
        result.point_add_assign(scaled);
        // Q_c * C
        scaled = vk.gate_setup_commitments[2].point_mul(proof.state_polys_openings_at_z[2]);
        result.point_add_assign(scaled);
        // Q_d * D
        scaled = vk.gate_setup_commitments[3].point_mul(proof.state_polys_openings_at_z[3]);
        result.point_add_assign(scaled);
        // Q_m* A*B or Q_ab*A*B
        PairingsBn254.Fr memory t = proof.state_polys_openings_at_z[0].copy();
        t.mul_assign(proof.state_polys_openings_at_z[1]);
        scaled = vk.gate_setup_commitments[4].point_mul(t);
        result.point_add_assign(scaled);
        // Q_AC* A*C
        t = proof.state_polys_openings_at_z[0].copy();
        t.mul_assign(proof.state_polys_openings_at_z[2]);
        scaled = vk.gate_setup_commitments[5].point_mul(t);
        result.point_add_assign(scaled);
        // Q_const
        result.point_add_assign(vk.gate_setup_commitments[6]);
        // Q_dNext * D_next
        scaled = vk.gate_setup_commitments[7].point_mul(proof.state_polys_openings_at_z_omega[0]);
        result.point_add_assign(scaled);        
        result.point_mul_assign(proof.gate_selectors_openings_at_z[0]);        

        PairingsBn254.G1Point memory rescue_custom_gate_linearization_contrib = rescue_custom_gate_linearization_contribution(vk, proof, state);
        result.point_add_assign(rescue_custom_gate_linearization_contrib);
        require(vk.non_residues.length == STATE_WIDTH-1);

        PairingsBn254.Fr memory one = PairingsBn254.new_fr(1);
        PairingsBn254.Fr memory factor = state.alpha_values[4].copy();
        for(uint256 i = 0;  i < proof.state_polys_openings_at_z.length; i++){
            t = state.z.copy();
            if(i == 0){
                t.mul_assign(one);
            }else{
                t.mul_assign(vk.non_residues[i-1]); // TODO add one into non-residues during codegen?
            }
            t.mul_assign(state.beta);
            t.add_assign(state.gamma);
            t.add_assign(proof.state_polys_openings_at_z[i]);

            factor.mul_assign(t);
        }

        scaled = proof.copy_permutation_grand_product_commitment.point_mul(factor);
        result.point_add_assign(scaled);
        
        // - (a(z) + beta*perm_a + gamma)*()*()*z(z*omega) * beta * perm_d(X)
        factor = state.alpha_values[4].copy();
        factor.mul_assign(state.beta);
        factor.mul_assign(proof.copy_permutation_grand_product_opening_at_z_omega);
        for(uint256 i = 0;  i < STATE_WIDTH-1; i++){
            t = proof.copy_permutation_polys_openings_at_z[i].copy();
            t.mul_assign(state.beta);
            t.add_assign(state.gamma);
            t.add_assign(proof.state_polys_openings_at_z[i]);

            factor.mul_assign(t);
        }
        scaled = vk.permutation_commitments[3].point_mul(factor);
        result.point_sub_assign(scaled);
        
        // + L_0(z) * Z(x)
        // TODO
        state.l_0_at_z = evaluate_lagrange_poly_out_of_domain(0, vk.domain_size, vk.omega, state.z);
        require(state.l_0_at_z.value != 0);
        factor = state.l_0_at_z.copy();
        factor.mul_assign(state.alpha_values[4 + 1]);
        scaled = proof.copy_permutation_grand_product_commitment.point_mul(factor);
        result.point_add_assign(scaled);

    }
    function rescue_custom_gate_linearization_contribution(VerificationKey memory vk, Proof memory proof, PartialVerifierState memory state) public view returns(PairingsBn254.G1Point memory result){        
        PairingsBn254.Fr memory t;
        PairingsBn254.Fr memory intermediate_result;

        // a^2 - b = 0
        t = proof.state_polys_openings_at_z[0].copy();
        t.mul_assign(t);
        t.sub_assign(proof.state_polys_openings_at_z[1]);
        // t.mul_assign(challenge1);
        t.mul_assign(state.alpha_values[1]);
        intermediate_result.add_assign(t);

        // b^2 - c = 0
        t = proof.state_polys_openings_at_z[1].copy();
        t.mul_assign(t);
        t.sub_assign(proof.state_polys_openings_at_z[2]);
        t.mul_assign(state.alpha_values[1+1]);
        intermediate_result.add_assign(t);

        // c*a - d = 0;
        t = proof.state_polys_openings_at_z[2].copy();
        t.mul_assign(proof.state_polys_openings_at_z[0]);
        t.sub_assign(proof.state_polys_openings_at_z[3]);
        t.mul_assign(state.alpha_values[1+2]);
        intermediate_result.add_assign(t);        

        result = vk.gate_selectors_commitments[1].point_mul(intermediate_result);        
    }
    struct Queries {
        PairingsBn254.G1Point[10] commitments_at_z;
        PairingsBn254.Fr[10] values_at_z;

        PairingsBn254.G1Point[3] commitments_at_z_omega;
        PairingsBn254.Fr[3] values_at_z_omega;
    }

    function prepare_queries(
        VerificationKey memory vk, 
        Proof memory proof, 
        PartialVerifierState memory state
    ) public view returns(Queries memory queries){
        // we set first two items in calee side so start idx from 2
        uint256 idx = 2;        
        for(uint256 i = 0; i<STATE_WIDTH; i++){
            queries.commitments_at_z[idx] = proof.state_polys_commitments[i];
            queries.values_at_z[idx] = proof.state_polys_openings_at_z[i];
            idx+=1;
        }
        require(proof.gate_selectors_openings_at_z.length == 1);
        queries.commitments_at_z[idx] = vk.gate_selectors_commitments[0];
        queries.values_at_z[idx] = proof.gate_selectors_openings_at_z[0];
        idx+=1;
        for(uint256 i = 0; i<STATE_WIDTH-1; i++){
            queries.commitments_at_z[idx] = vk.permutation_commitments[i];
            queries.values_at_z[idx] = proof.copy_permutation_polys_openings_at_z[i];
            idx+=1;
        }

        queries.commitments_at_z_omega[0] = proof.copy_permutation_grand_product_commitment;
        queries.commitments_at_z_omega[1] = proof.state_polys_commitments[STATE_WIDTH-1];

        queries.values_at_z_omega[0] = proof.copy_permutation_grand_product_opening_at_z_omega;
        queries.values_at_z_omega[1] = proof.state_polys_openings_at_z_omega[0];

    }
    
    function final_pairing(
        // VerificationKey memory vk, 
        PairingsBn254.G2Point[NUM_G2_ELS] memory g2_elements,
        Proof memory proof, 
        PartialVerifierState memory state,
        PairingsBn254.G1Point memory aggregated_commitment_at_z,
        PairingsBn254.G1Point memory aggregated_commitment_at_z_omega,
        PairingsBn254.Fr memory aggregated_opening_at_z,
        PairingsBn254.Fr memory aggregated_opening_at_z_omega
        ) internal view returns(bool){

        // q(x) = f(x) - f(z) / (x - z)
        // q(x) * (x-z)  = f(x) - f(z)

        // f(x)
        PairingsBn254.G1Point memory  pair_with_generator = aggregated_commitment_at_z.copy_g1();
        aggregated_commitment_at_z_omega.point_mul_assign(state.u);
        pair_with_generator.point_add_assign(aggregated_commitment_at_z_omega);

        // - f(z)*g
        PairingsBn254.Fr memory  aggregated_value = aggregated_opening_at_z_omega.copy();
        aggregated_value.mul_assign(state.u);
        aggregated_value.add_assign(aggregated_opening_at_z);
        PairingsBn254.G1Point memory  tp = PairingsBn254.P1().point_mul(aggregated_value);
        pair_with_generator.point_sub_assign(tp);

        // +z * q(x)
        tp = proof.opening_proof_at_z.point_mul(state.z);
        PairingsBn254.Fr memory t = state.z_omega.copy();
        t.mul_assign(state.u);
        PairingsBn254.G1Point memory t1 = proof.opening_proof_at_z_omega.point_mul(t);
        tp.point_add_assign(t1);
        pair_with_generator.point_add_assign(tp);

        // rhs
        PairingsBn254.G1Point memory pair_with_x = proof.opening_proof_at_z_omega.point_mul(state.u);
        pair_with_x.point_add_assign(proof.opening_proof_at_z);
        pair_with_x.negate();
        // Pairing precompile expects points to be in a `i*x[1] + x[0]` form instead of `x[0] + i*x[1]`
        // so we handle it in code generation step
        PairingsBn254.G2Point memory first_g2 = g2_elements[0];
        PairingsBn254.G2Point memory second_g2 = g2_elements[1];
        PairingsBn254.G2Point memory gen2 = PairingsBn254.P2();
                
        bool valid = PairingsBn254.pairingProd2(pair_with_generator, first_g2, pair_with_x, second_g2);         
        return valid;               
    }
}


contract Verifier is Plonk4VerifierWithAccessToDNext{

    function get_verification_key() internal pure returns(VerificationKey memory vk){
        vk.num_inputs = 0;
        vk.domain_size = 512;
        vk.omega = PairingsBn254.new_fr(0x0dd30b9ad8c173555d2a33029bc807ac165b61281e9054a173af7ff4e4fc88fc);
        // coefficients
        vk.gate_setup_commitments[0] = PairingsBn254.new_g1(0x1fb1a6158a51c9cfd190c2a96846d99c60a9982be8e5c9a9e12ff7ba0b54f2a1,0x24a0119fdfab12365fa9295aa905565c007c15b2e6fe2872882b3cce393c091c);
        vk.gate_setup_commitments[1] = PairingsBn254.new_g1(0x2a03abd51f0145565abbe3f9c031f9b91b8eec5926c85a3a8db243a46cd3d952,0x29ede32d47094c2c1d3c75641701e2476d57f260250cb6f2f67bb7ff0c5ec3af);
        vk.gate_setup_commitments[2] = PairingsBn254.new_g1(0x0b692a950e9461f578f15cc0d2e608b250ae56bbd498d771eb3d656c77747844,0x147fc64f1038b2e19d5e33972416ce2e6dcce1c92048db5d02781c2dca63be33);
        vk.gate_setup_commitments[3] = PairingsBn254.new_g1(0x14e7adde62def734986886a6412aa519a5a70e5a36efd9a3ef92919f6a5c41d4,0x093f801bb05783a8855421826428575af5d154ded16eb116c0cf842ef4fa01fa);
        vk.gate_setup_commitments[4] = PairingsBn254.new_g1(0x2a0fc50778a60980d36b1c47bad8217093cc06adbe66795d375628f59bf97bbc,0x06241e9d2bb1b11c986d8fbfde8bc92cf3783a12ced73071c06b146fc99a97bf);
        vk.gate_setup_commitments[5] = PairingsBn254.new_g1(0x1f5d769347d3cba40d02f91178d588735a9c457d52d354fd7bc95d7527dce43f,0x176a340d5fea07544a37edffc43a6c9b7c944604f63fa2f7bfe5d427e41942e7);
        vk.gate_setup_commitments[6] = PairingsBn254.new_g1(0x118267d45429c9ee512b7e5b53e16e80c9b5e49c24a319fbf165fc9a641861d5,0x195d2f076801781936637d810e0627ac28fa5ba9c302f274f7eb073c57d5fea9);
        vk.gate_setup_commitments[7] = PairingsBn254.new_g1(0x0000000000000000000000000000000000000000000000000000000000000000,0x0000000000000000000000000000000000000000000000000000000000000001);
        // gate selectors
        vk.gate_selectors_commitments[0] = PairingsBn254.new_g1(0x1426778b2e8c3e26b2f5b535905c5478903bf1e383587bdaa01e53aed87fd9e3,0x08572ee7826a7bb7fd37e6811ed99f496a75eae44be75f1f22e87ff2560c1b19);
        vk.gate_selectors_commitments[1] = PairingsBn254.new_g1(0x14c24f853bf75e6e96352d9980dbe40b09134118824d622d4818415711c3157d,0x28f20057ee04a164cf7e062b7ea65adf106af714c5a615c4ccff383709c275d8);
        // permutation
        vk.permutation_commitments[0] = PairingsBn254.new_g1(0x0f258c86c04a30b55ae13ed70cbcdd0ae4394bfa554f9d0d4d4cb3f0921546e8,0x25cdbcfa477ad6d74b3ea1678a1691e126aedd9c5e5de3139bcc92e3876d8bd0);
        vk.permutation_commitments[1] = PairingsBn254.new_g1(0x01b8836eae18b732c2d79ce8bf2575a894764132391e44c038bc9d00911f5d5f,0x117ff45a60d7f41c238cf7e3cee8981346912267c77854b023ce7c464d157c6d);
        vk.permutation_commitments[2] = PairingsBn254.new_g1(0x2242faa0b9f1b72856fa7337f2eaa2747510c71c12aeb3a0acf00a1e1115e29e,0x29c94a75a018f556335fb36eb2708165f390945aaa9dc49553ceea5915c87972);
        vk.permutation_commitments[3] = PairingsBn254.new_g1(0x0268c1d0b316faec267ffd34beffac936d08a1716b60d9ac7003fc9cad3aff06,0x194eaf0a8c6b63948d1ad8dae6c180288b4525bee0baafb03b785bae7eaa0159);
        // non residues
        vk.non_residues[0] = PairingsBn254.new_fr(0x0000000000000000000000000000000000000000000000000000000000000005);
        vk.non_residues[1] = PairingsBn254.new_fr(0x0000000000000000000000000000000000000000000000000000000000000007);
        vk.non_residues[2] = PairingsBn254.new_fr(0x000000000000000000000000000000000000000000000000000000000000000a);
        
        // g2 elements
        vk.g2_elements[0] = PairingsBn254.new_g2([0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2,0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed],[0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b,0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa]);
        vk.g2_elements[1] = PairingsBn254.new_g2([0x12740934ba9615b77b6a49b06fcce83ce90d67b1d0e2a530069e3a7306569a91,0x116da8c89a0d090f3d8644ada33a5f1c8013ba7204aeca62d66d931b99afe6e7],[0x25222d9816e5f86b4a7dedd00d04acc5c979c18bd22b834ea8c6d07c0ba441db,0x076441042e77b6309644b56251f059cf14befc72ac8a6157d30924e58dc4c172]);
    }

    function deserialize_proof(
        uint256[] memory public_inputs, 
        uint256[] memory serialized_proof
    ) internal view returns(Proof memory proof) {
        // require(serialized_proof.length == 44); TODO
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
        for (uint256 i = 0; i < proof.quotient_poly_parts_commitments.length; i++) {
            proof.quotient_poly_parts_commitments[i] = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
            );
            j += 2;
        }

        for (uint256 i = 0; i < proof.state_polys_openings_at_z.length; i++) {
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


library PairingsBn254 {
    uint256 constant q_mod = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
    uint256 constant r_mod = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
    uint256 constant bn254_b_coeff = 3;

    struct G1Point {
        uint256 X;
        uint256 Y;
    } 
    
    struct Fr {
        uint256 value;
    }
    
    function new_fr(uint256 fr) internal pure returns (Fr memory) {
        require(fr < r_mod);
        return Fr({value: fr});
    }
    
    function copy(Fr memory self) internal pure returns (Fr memory n) {
        n.value = self.value;
    }
    
    function assign(Fr memory self, Fr memory other) internal pure {
        self.value = other.value;
    }
    
    function inverse(Fr memory fr) internal view returns (Fr memory) {
        require(fr.value != 0);
        return pow(fr, r_mod-2);
    }
    
    function add_assign(Fr memory self, Fr memory other) internal pure {
        self.value = addmod(self.value, other.value, r_mod);
    }
    
    function sub_assign(Fr memory self, Fr memory other) internal pure {
        self.value = addmod(self.value, r_mod - other.value, r_mod);
    }
    
    function mul_assign(Fr memory self, Fr memory other) internal pure {
        self.value = mulmod(self.value, other.value, r_mod);
    }
    
    function pow(Fr memory self, uint256 power) internal view returns (Fr memory) {
        uint256[6] memory input = [32, 32, 32, self.value, power, r_mod];
        uint256[1] memory result;
        bool success;
        assembly {
            success := staticcall(gas(), 0x05, input, 0xc0, result, 0x20)
        }
        require(success);
        return Fr({value: result[0]});
    }
    
    // Encoding of field elements is: X[0] * z + X[1]
    struct G2Point {
        uint[2] X;
        uint[2] Y;
    }

    function P1() internal pure returns (G1Point memory) {
        return G1Point(1, 2);
    }
    
    function new_g1(uint256 x, uint256 y) internal pure returns (G1Point memory) {
        return G1Point(x, y);
    }

    // function new_g1_checked(uint256 x, uint256 y) internal pure returns (G1Point memory) {
    function new_g1_checked(uint256 x, uint256 y) internal pure returns (G1Point memory) {
        if (x == 0 && y == 0) {
            // point of infinity is (0,0)
            return G1Point(x, y);
        }
        
        // check encoding
        require(x < q_mod, "x axis isn't valid");
        require(y < q_mod, "y axis isn't valid");
        // check on curve
        uint256 lhs = mulmod(y, y, q_mod); // y^2
        
        uint256 rhs = mulmod(x, x, q_mod); // x^2
        rhs = mulmod(rhs, x, q_mod); // x^3        
        rhs = addmod(rhs, bn254_b_coeff, q_mod); // x^3 + b
        require(lhs == rhs, "is not on curve");

        return G1Point(x, y);
    }
    
    function new_g2(uint256[2] memory x, uint256[2] memory y) internal pure returns (G2Point memory) {
        return G2Point(x, y);
    }
    
    function copy_g1(G1Point memory self) internal pure returns (G1Point memory result) {
        result.X = self.X;
        result.Y = self.Y;
    }

    function P2() internal pure returns (G2Point memory) {
        // for some reason ethereum expects to have c1*v + c0 form
        
        return G2Point(
            [0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2,
                0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed],
            [0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b,
                0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa]
        );
    }

    function negate(G1Point memory self) internal pure {
        // The prime q in the base field F_q for G1
        if (self.Y == 0) {
            require(self.X == 0);
            return;
        }

        self.Y = q_mod - self.Y;
    }

    function point_add(G1Point memory p1, G1Point memory p2)
        internal view returns (G1Point memory r)
    {
        point_add_into_dest(p1, p2, r);
        return r;
    }
    
    function point_add_assign(G1Point memory p1, G1Point memory p2)
        internal view
    {
        point_add_into_dest(p1, p2, p1);
    }

    function point_add_into_dest(G1Point memory p1, G1Point memory p2, G1Point memory dest)
        internal view
    {
        if (p2.X == 0 && p2.Y == 0) {
            // we add zero, nothing happens
            dest.X = p1.X;
            dest.Y = p1.Y;
            return;
        } else if (p1.X == 0 && p1.Y == 0) {
            // we add into zero, and we add non-zero point
            dest.X = p2.X;
            dest.Y = p2.Y;
            return;
        } else {
            uint256[4] memory input;

            input[0] = p1.X;
            input[1] = p1.Y;
            input[2] = p2.X;
            input[3] = p2.Y;

            bool success = false;
            assembly {
                success := staticcall(gas(), 6, input, 0x80, dest, 0x40)
            }
            require(success);
        }
    }
    
    function point_sub_assign(G1Point memory p1, G1Point memory p2)
        internal view
    {
        point_sub_into_dest(p1, p2, p1);
    }

    function point_sub_into_dest(G1Point memory p1, G1Point memory p2, G1Point memory dest)
        internal view
    {
        if (p2.X == 0 && p2.Y == 0) {
            // we subtracted zero, nothing happens
            dest.X = p1.X;
            dest.Y = p1.Y;
            return;
        } else if (p1.X == 0 && p1.Y == 0) {
            // we subtract from zero, and we subtract non-zero point
            dest.X = p2.X;
            dest.Y = q_mod - p2.Y;
            return;
        } else {
            uint256[4] memory input;

            input[0] = p1.X;
            input[1] = p1.Y;
            input[2] = p2.X;
            input[3] = q_mod - p2.Y;

            bool success = false;
            assembly {
                success := staticcall(gas(), 6, input, 0x80, dest, 0x40)
            }
            require(success);
        }
    }

    function point_mul(G1Point memory p, Fr memory s)
        internal view returns (G1Point memory r)
    {
        // https://eips.ethereum.org/EIPS/eip-197
        // Elliptic curve points are encoded as a Jacobian pair (X, Y) where the point at infinity is encoded as (0, 0)
        // TODO
        if(p.X == 0 && p.Y == 1){
            p.Y = 0;
        }
        point_mul_into_dest(p, s, r);
        return r;
    }
    
    function point_mul_assign(G1Point memory p, Fr memory s)
        internal view
    {
        point_mul_into_dest(p, s, p);
    }

    function point_mul_into_dest(G1Point memory p, Fr memory s, G1Point memory dest)
        internal view
    {
        uint[3] memory input;   
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s.value;
        bool success;
        assembly {
            success := staticcall(gas(), 7, input, 0x60, dest, 0x40)
        }
        require(success);
    }
    
    function pairing(G1Point[] memory p1, G2Point[] memory p2)
        internal view returns (bool)
    {
        require(p1.length == p2.length);
        uint elements = p1.length;
        uint inputSize = elements * 6;
        uint[] memory input = new uint[](inputSize);
        for (uint i = 0; i < elements; i++)
        {
            input[i * 6 + 0] = p1[i].X;
            input[i * 6 + 1] = p1[i].Y;
            input[i * 6 + 2] = p2[i].X[0];
            input[i * 6 + 3] = p2[i].X[1];
            input[i * 6 + 4] = p2[i].Y[0];
            input[i * 6 + 5] = p2[i].Y[1];
        }
        uint[1] memory out;
        bool success;
        assembly {
            success := staticcall(gas(), 8, add(input, 0x20), mul(inputSize, 0x20), out, 0x20)
        }
        require(success);
        return out[0] != 0;
    }

    /// Convenience method for a pairing check for two pairs.
    function pairingProd2(G1Point memory a1, G2Point memory a2, G1Point memory b1, G2Point memory b2)
        internal view returns (bool)
    {
        G1Point[] memory p1 = new G1Point[](2);
        G2Point[] memory p2 = new G2Point[](2);
        p1[0] = a1;
        p1[1] = b1;
        p2[0] = a2;
        p2[1] = b2;
        return pairing(p1, p2);
    }
}


library TranscriptLibrary {
    // flip                    0xe000000000000000000000000000000000000000000000000000000000000000;
    uint256 constant FR_MASK = 0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff;

    uint32 constant DST_0 = 0;
    uint32 constant DST_1 = 1;
    uint32 constant DST_CHALLENGE = 2;
    
    struct Transcript {
        bytes32 state_0;
        bytes32 state_1;
        uint32 challenge_counter;
    }

    function new_transcript() internal pure returns (Transcript memory t) {
        t.state_0 = bytes32(0);
        t.state_1 = bytes32(0);
        t.challenge_counter = 0;
    }

    function update_with_u256(Transcript memory self, uint256 value) internal pure {
        bytes32 old_state_0 = self.state_0;
        self.state_0 = keccak256(abi.encodePacked(DST_0, old_state_0, self.state_1, value));
        self.state_1 = keccak256(abi.encodePacked(DST_1, old_state_0, self.state_1, value));
    }
    
    function update_with_fr(Transcript memory self, PairingsBn254.Fr memory value) internal pure {
        update_with_u256(self, value.value);
    }
    
    function update_with_g1(Transcript memory self, PairingsBn254.G1Point memory p) internal pure {
        update_with_u256(self, p.X);
        update_with_u256(self, p.Y);
    }
    
    function get_challenge(Transcript memory self) internal pure returns(PairingsBn254.Fr memory challenge) {
        bytes32 query = keccak256(abi.encodePacked(DST_CHALLENGE, self.state_0, self.state_1, self.challenge_counter));
        self.challenge_counter += 1;
        challenge = PairingsBn254.Fr({value: uint256(query) & FR_MASK});
    }
}
