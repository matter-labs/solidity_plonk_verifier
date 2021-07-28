pragma solidity >0.7;

import "./Pairing.sol";
import "./Transcript.sol";
import "./VerificationKey.sol";

uint256 constant SW = 4; // TODO

contract Plonk4VerifierWithAccessToDNext {
    using PairingsBn254 for PairingsBn254.G1Point;
    using PairingsBn254 for PairingsBn254.G2Point;
    using PairingsBn254 for PairingsBn254.Fr;
    
    using TranscriptLibrary for TranscriptLibrary.Transcript;
    
    struct Proof {
        uint256[] input_values;
        // commitments
        PairingsBn254.G1Point[SW] state_polys_commitments;
        PairingsBn254.G1Point copy_permutation_grand_product_commitment;
        PairingsBn254.G1Point[SW] quotient_poly_parts_commitments;
        PairingsBn254.G1Point lookup_s_poly_commitment;
        PairingsBn254.G1Point lookup_grand_product_commitment;
        // openings
        PairingsBn254.Fr[SW] state_polys_openings_at_z;
        PairingsBn254.Fr[1] state_polys_openings_at_z_omega; // TODO: only D_next
        PairingsBn254.Fr[] gate_setup_openings_at_z;
        PairingsBn254.Fr[] gate_selectors_openings_at_z;
        PairingsBn254.Fr[] copy_permutation_polys_openings_at_z;
        PairingsBn254.Fr copy_permutation_grand_product_opening_at_z_omega;
        PairingsBn254.Fr quotient_poly_opening_at_z;
        PairingsBn254.Fr linearization_poly_opening_at_z;
        PairingsBn254.Fr[SW-1] permutation_polynomials_at_z;
        // lookup openings
        PairingsBn254.Fr lookup_s_poly_opening_at_z_omega;
        PairingsBn254.Fr lookup_grand_product_opening_at_z_omega;
        PairingsBn254.Fr lookup_t_poly_opening_at_z;
        PairingsBn254.Fr lookup_t_poly_opening_at_z_omega;
        PairingsBn254.Fr lookup_selector_poly_opening_at_z;
        PairingsBn254.Fr lookup_table_type_poly_opening_at_z;

    
        PairingsBn254.G1Point opening_proof_at_z;
        PairingsBn254.G1Point opening_proof_at_z_omega;
    }
    
    struct PartialVerifierState {
        PairingsBn254.Fr eta;
        PairingsBn254.Fr alpha;
        PairingsBn254.Fr alpha4;
        PairingsBn254.Fr alpha5;
        PairingsBn254.Fr alpha6;
        PairingsBn254.Fr alpha7;
        PairingsBn254.Fr alpha8;
        PairingsBn254.Fr beta;
        PairingsBn254.Fr gamma;
        PairingsBn254.Fr beta_lookup;
        PairingsBn254.Fr gamma_lookup;
        PairingsBn254.Fr beta_plus_one;
        PairingsBn254.Fr beta_gamma;
        PairingsBn254.Fr v;
        PairingsBn254.Fr u;
        PairingsBn254.Fr z;
        PairingsBn254.Fr[] cached_lagrange_evals;
        PairingsBn254.Fr z_minus_last_omega;

    }

    
    function evaluate_l0_at_point(
        uint256 domain_size, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr memory num) {
        PairingsBn254.Fr memory one = PairingsBn254.new_fr(1);

        PairingsBn254.Fr memory size_fe = PairingsBn254.new_fr(domain_size);
        PairingsBn254.Fr memory den = at;
        den.sub_assign(one);
        den.mul_assign(size_fe);

        den = den.inverse(); // TODO: check has inverse

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
    
    function batch_evaluate_lagrange_poly_out_of_domain(
        uint256[] memory poly_nums, 
        uint256 domain_size, 
        PairingsBn254.Fr memory omega, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr[] memory res) {
        PairingsBn254.Fr memory one = PairingsBn254.new_fr(1);
        PairingsBn254.Fr memory tmp_1 = PairingsBn254.new_fr(0);
        PairingsBn254.Fr memory tmp_2 = PairingsBn254.new_fr(domain_size);
        PairingsBn254.Fr memory vanishing_at_z = at.pow(domain_size);
        vanishing_at_z.sub_assign(one);
        // we can not have random point z be in domain
        require(vanishing_at_z.value != 0);
        PairingsBn254.Fr[] memory nums = new PairingsBn254.Fr[](poly_nums.length);
        PairingsBn254.Fr[] memory dens = new PairingsBn254.Fr[](poly_nums.length);
        // numerators in a form omega^i * (z^n - 1)
        // denoms in a form (z - omega^i) * N
        for (uint i = 0; i < poly_nums.length; i++) {
            tmp_1 = omega.pow(poly_nums[i]); // power of omega
            nums[i].assign(vanishing_at_z);
            nums[i].mul_assign(tmp_1);
            
            dens[i].assign(at); // (X - omega^i) * N
            dens[i].sub_assign(tmp_1); 
            dens[i].mul_assign(tmp_2); // mul by domain size
        }
        
        PairingsBn254.Fr[] memory partial_products = new PairingsBn254.Fr[](poly_nums.length);
        partial_products[0].assign(PairingsBn254.new_fr(1));
        for (uint i = 1; i < dens.length; i++) {
            partial_products[i].assign(dens[i-1]);
            partial_products[i].mul_assign(partial_products[i-1]);
        }
    
        tmp_2.assign(partial_products[partial_products.length - 1]);
        tmp_2.mul_assign(dens[dens.length - 1]);
        tmp_2 = tmp_2.inverse(); // tmp_2 contains a^-1 * b^-1 (with! the last one)
        
        for (uint i = dens.length - 1; i < dens.length; i--) {
            tmp_1.assign(tmp_2); // all inversed
            tmp_1.mul_assign(partial_products[i]); // clear lowest terms
            tmp_2.mul_assign(dens[i]);
            dens[i].assign(tmp_1);
        }
        
        for (uint i = 0; i < nums.length; i++) {
            nums[i].mul_assign(dens[i]);
        }

        return nums;
    }
    
    function evaluate_vanishing(
        uint256 domain_size, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr memory res) {
        res = at.pow(domain_size);
        res.sub_assign(PairingsBn254.new_fr(1));
    }
    
    function verify(Proof memory proof, VerificationKey memory vk) internal view returns (bool) {
        PartialVerifierState memory state;

        TranscriptLibrary.Transcript memory transcript = TranscriptLibrary.new_transcript();

        for(uint256 i =0; i < vk.num_inputs; i++){
            transcript.update_with_u256(proof.input_values[i]);
        }

        for(uint256 i =0; i < SW; i++){
            transcript.update_with_g1(proof.state_polys_commitments[i]);
        }

        bool has_lookup = vk.total_lookup_entries_length > 0;

        if(has_lookup){
            state.eta = transcript.get_challenge();
            transcript.update_with_g1(proof.lookup_s_poly_commitment);
        }

        transcript.update_with_g1(proof.copy_permutation_grand_product_commitment);

        if(has_lookup){
            state.beta_lookup = transcript.get_challenge();
            state.gamma_lookup = transcript.get_challenge();

            transcript.update_with_g1(proof.lookup_grand_product_commitment);
        }

        state.alpha = transcript.get_challenge();

        for(uint256 i =0; i < proof.quotient_poly_parts_commitments.length; i++){
            transcript.update_with_g1(proof.quotient_poly_parts_commitments[i]);
        }   

        state.z = transcript.get_challenge();

        transcript.update_with_fr(proof.quotient_poly_opening_at_z);

        for(uint256 i =0; i < proof.state_polys_openings_at_z.length; i++){
            transcript.update_with_fr(proof.state_polys_openings_at_z[i]);
        }          

        for(uint256 i =0; i < proof.copy_permutation_polys_openings_at_z.length; i++){
            transcript.update_with_fr(proof.copy_permutation_polys_openings_at_z[i]);
        }          

        PairingsBn254.Fr memory z_omega = state.z;
        z_omega.mul_assign(vk.omega);

        transcript.update_with_fr(proof.copy_permutation_grand_product_opening_at_z_omega);

        if(has_lookup){
            transcript.update_with_fr(proof.lookup_t_poly_opening_at_z);
            transcript.update_with_fr(proof.lookup_selector_poly_opening_at_z);
            transcript.update_with_fr(proof.lookup_table_type_poly_opening_at_z);
            transcript.update_with_fr(proof.lookup_s_poly_opening_at_z_omega);
            transcript.update_with_fr(proof.lookup_grand_product_opening_at_z_omega);
            transcript.update_with_fr(proof.lookup_t_poly_opening_at_z_omega);
        }
        transcript.update_with_fr(proof.linearization_poly_opening_at_z);

        state.v = transcript.get_challenge();
        
        transcript.update_with_g1(proof.opening_proof_at_z);
        transcript.update_with_g1(proof.opening_proof_at_z_omega);

        state.u = transcript.get_challenge();

        
        uint256[] memory lagrange_poly_numbers = new uint256[](vk.num_inputs);
        for (uint256 i = 0; i < lagrange_poly_numbers.length; i++) {
            lagrange_poly_numbers[i] = i;
        }
        
        state.cached_lagrange_evals = batch_evaluate_lagrange_poly_out_of_domain(
            lagrange_poly_numbers,
            vk.domain_size, 
            vk.omega, 
            state.z
        );


        require(vk.num_inputs > 0);
        
        PairingsBn254.Fr memory tmp;
        PairingsBn254.Fr memory inputs_term;
        for(uint256 i =0; i < vk.num_inputs; i++){
            tmp = state.cached_lagrange_evals[i];
            tmp.mul_assign(PairingsBn254.new_fr(proof.input_values[i]));

            inputs_term.add_assign(tmp);
        }          
        inputs_term.mul_assign(PairingsBn254.new_fr(1));

        PairingsBn254.Fr memory t_num_on_full_domain = proof.gate_selectors_openings_at_z[0];
        t_num_on_full_domain.mul_assign(inputs_term);
        t_num_on_full_domain.add_assign(proof.linearization_poly_opening_at_z);


        // alpha values for permutation argument
        state.alpha4 = state.alpha.pow(4);
        state.alpha5 = state.alpha4;
        state.alpha5.mul_assign(state.alpha);

        PairingsBn254.Fr memory factor = state.alpha4;
        factor.mul_assign(proof.copy_permutation_grand_product_opening_at_z_omega);

        require(proof.copy_permutation_polys_openings_at_z.length == 3);

        PairingsBn254.Fr memory t;  // TMP;
        // - alpha_0 * (a + perm(z) * beta + gamma)*()*(d + gamma) * z(z*omega)
        for(uint256 i = 0; i < proof.copy_permutation_polys_openings_at_z.length; i++){
            t = proof.copy_permutation_polys_openings_at_z[i];
            t.mul_assign(state.beta);
            t.add_assign(proof.state_polys_openings_at_z[i]);
            t.add_assign(state.gamma);

            factor.mul_assign(t);
        }

        require(SW == 4);

        t = proof.state_polys_openings_at_z[3];
        t.add_assign(state.gamma);
        factor.mul_assign(t);

        t_num_on_full_domain.sub_assign(factor);

        // - L_0(z) * alpha_1
        // TODO check domain is power of two
        PairingsBn254.Fr memory l_0_at_z = evaluate_l0_at_point(vk.domain_size, state.z);
        l_0_at_z.mul_assign(state.alpha5);

        t_num_on_full_domain.sub_assign(l_0_at_z);

        PairingsBn254.Fr memory zero = PairingsBn254.new_fr(0);
        PairingsBn254.Fr memory one = PairingsBn254.new_fr(1);
        if(has_lookup){
            state.alpha6 = state.alpha.pow(6);
            state.alpha7 = state.alpha6;
            state.alpha7.mul_assign(state.alpha);
            state.alpha8 = state.alpha7;
            state.alpha8.mul_assign(state.alpha);

            state.beta_plus_one = state.beta_lookup;
            state.beta_plus_one.add_assign(one);

            state.beta_gamma = state.beta_plus_one;
            state.beta_gamma.mul_assign(state.gamma_lookup);

            // (s'*beta + gamma)*(zw')*alpha
            t = proof.lookup_s_poly_opening_at_z_omega;
            t.mul_assign(state.beta_lookup);
            t.add_assign(state.beta_gamma);
            t.mul_assign(proof.lookup_grand_product_opening_at_z_omega);
            t.mul_assign(state.alpha6);

            // (z - omega^{n-1}) for this part
            PairingsBn254.Fr memory last_omega = vk.omega.pow(vk.domain_size -1);
            state.z_minus_last_omega = state.z;
            state.z_minus_last_omega.sub_assign(last_omega);
            t.mul_assign(state.z_minus_last_omega);

            t_num_on_full_domain.add_assign(t);

            // - alpha_1 * L_{0}(z)
            t = l_0_at_z;
            t.mul_assign(state.alpha7);        
            t_num_on_full_domain.sub_assign(t);

            // - alpha_2 * beta_gamma_powered L_{n-1}(z)
            PairingsBn254.Fr memory beta_gamma_powered = state.beta_gamma.pow(vk.domain_size-1);
            // TODO
            PairingsBn254.Fr memory l_n_minus_one_at_z = evaluate_lagrange_poly_out_of_domain(0, vk.domain_size -1, vk.omega, state.z);
            l_n_minus_one_at_z.mul_assign(beta_gamma_powered);
            l_n_minus_one_at_z.mul_assign(state.alpha8);
            
            t_num_on_full_domain.sub_assign(l_n_minus_one_at_z);

            // sanity check
            PairingsBn254.Fr memory lhs = proof.quotient_poly_opening_at_z;
            lhs.mul_assign(evaluate_vanishing(vk.domain_size, state.z));
            lhs.sub_assign(t_num_on_full_domain);
            require(lhs.value == 0);
        }


        require(vk.gate_setup_commitments.length == 7);
        require(proof.state_polys_openings_at_z.length == 4);
        require(proof.state_polys_openings_at_z_omega.length == 1);
        // rescue custom gate doesn't contribute to linearization, 
        // only gate selectors participate into verification
        // qMain*(Q_a * A + Q_b * B + Q_c * C + Q_d * D + Q_m * A*B + Q_const + Q_dNext * D_next)
        PairingsBn254.G1Point memory r = PairingsBn254.new_g1(0, 0);
        // Q_a * A        
        PairingsBn254.G1Point memory scaled = vk.gate_setup_commitments[0].point_mul(proof.state_polys_openings_at_z[0]);
        r.point_add_assign(scaled);

        // Q_b * B
        scaled = vk.gate_setup_commitments[1].point_mul(proof.state_polys_openings_at_z[1]);
        r.point_add_assign(scaled);

        // Q_c * C
        scaled = vk.gate_setup_commitments[2].point_mul(proof.state_polys_openings_at_z[2]);
        r.point_add_assign(scaled);

        // Q_d * D
        scaled = vk.gate_setup_commitments[3].point_mul(proof.state_polys_openings_at_z[3]);
        r.point_add_assign(scaled);

        // Q_m* A*B
        t = proof.state_polys_openings_at_z[0];
        t.mul_assign(proof.state_polys_openings_at_z[1]);
        scaled = vk.gate_setup_commitments[4].point_mul(t);
        r.point_add_assign(scaled);

        // Q_const
        r.point_add_assign(vk.gate_setup_commitments[5]);

        // Q_dNext * D_next
        scaled = vk.gate_setup_commitments[6].point_mul(proof.state_polys_openings_at_z_omega[0]);
        r.point_add_assign(scaled);

        r.point_mul_assign(proof.gate_selectors_openings_at_z[0]);

        // now proceed with rescue custom gate, it only involves with selectors
        // we need 3 challenges
        PairingsBn254.Fr memory challenge_1 = state.alpha;
        PairingsBn254.Fr memory challenge_2 = challenge_1;
        challenge_2.mul_assign(state.alpha);
        PairingsBn254.Fr memory challenge_3 = challenge_2;
        challenge_3.mul_assign(state.alpha);

        PairingsBn254.Fr memory result;

        // a^2 - b = 0
        t = proof.state_polys_openings_at_z[0];
        t.mul_assign(t);
        t.sub_assign(proof.state_polys_openings_at_z[1]);
        result.add_assign(t);

        // b^2 - c = 0
        t = proof.state_polys_openings_at_z[1];
        t.mul_assign(t);
        t.sub_assign(proof.state_polys_openings_at_z[2]);
        result.add_assign(t);

        // c*a - d = 0;
        t = proof.state_polys_openings_at_z[2];
        t.mul_assign(proof.state_polys_openings_at_z[0]);
        t.sub_assign(proof.state_polys_openings_at_z[3]);
        result.add_assign(t);

        scaled = vk.gate_selectors_commitments[1].point_mul(result);
        r.point_add_assign(scaled);

        require(vk.non_residues.length == SW-1);
        // reuse factor
        factor = state.alpha4;
        for(uint256 i = 0;  i < proof.state_polys_openings_at_z.length; i++){
            t = state.z;
            if(i == 0){
                t.mul_assign(one);
            }else{
                t.mul_assign(vk.non_residues[i]);
            }
            t.mul_assign(state.beta);
            t.mul_assign(state.gamma);
            t.mul_assign(proof.state_polys_openings_at_z[i]);

            factor.mul_assign(t);
        }

        scaled = proof.copy_permutation_grand_product_commitment.point_mul(factor);
        r.point_add_assign(scaled);

        // - (a(z) + beta*perm_a + gamma)*()*()*z(z*omega) * beta * perm_d(X)
        factor = state.alpha4;
        factor.mul_assign(state.beta);
        factor.mul_assign(proof.copy_permutation_grand_product_opening_at_z_omega);
        for(uint256 i = 0;  i < 3; i++){
            t = proof.copy_permutation_polys_openings_at_z[i];
            t.mul_assign(state.beta);
            t.add_assign(state.gamma);
            t.add_assign(proof.state_polys_openings_at_z[i]);

            factor.mul_assign(t);
        }
        scaled = vk.permutation_commitments[3].point_mul(factor);
        r.point_sub_assign(scaled);

        // + L_0(z) * Z(x)
        factor = l_0_at_z;
        factor.mul_assign(state.alpha5);
        scaled = proof.copy_permutation_grand_product_commitment.point_mul(factor);
        r.point_add_assign(scaled);

        if(has_lookup){
            // s(x) from the Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega)))
            factor = proof.lookup_grand_product_opening_at_z_omega;
            factor.mul_assign(state.alpha6);
            factor.mul_assign(state.z_minus_last_omega);
            
            scaled = proof.lookup_s_poly_commitment.point_mul(factor);
            r.point_add_assign(scaled);

            // Z(x) from - alpha_0 * Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) 
            // + alpha_1 * Z(x) * L_{0}(z) + alpha_2 * Z(x) * L_{n-1}(z)

            // accumulate coefficient
            factor = proof.lookup_t_poly_opening_at_z_omega;
            factor.mul_assign(state.beta_lookup);
            factor.add_assign(proof.lookup_t_poly_opening_at_z);
            factor.add_assign(state.beta_gamma);

            // (\gamma + f(x))
            PairingsBn254.Fr memory f_reconstructed; // TODO tmp
            PairingsBn254.Fr memory current = one;
            PairingsBn254.Fr memory tmp0;
            for(uint256 i=0; i< proof.state_polys_openings_at_z.length; i++){
                tmp0 = proof.state_polys_openings_at_z[i];
                tmp0.mul_assign(current);
                f_reconstructed.add_assign(tmp0);

                current.mul_assign(state.eta);
            }

            // add type of table
            t = proof.lookup_table_type_poly_opening_at_z;
            t.mul_assign(current);
            f_reconstructed.mul_assign(t);

            f_reconstructed.mul_assign(proof.lookup_selector_poly_opening_at_z);
            f_reconstructed.add_assign(state.gamma_lookup);

            // end of (\gamma + f(x)) part
            factor.mul_assign(f_reconstructed);
            factor.mul_assign(state.beta_plus_one);
            t = zero;
            t.sub_assign(factor);
            factor = t;
            factor.mul_assign(state.alpha6);

            // Multiply by (z - omega^{n-1})
            factor.mul_assign(state.z_minus_last_omega);

            // L_{0}(z) in front of Z(x)
            t = l_0_at_z;
            t.mul_assign(state.alpha7);
            factor.add_assign(t);

            scaled = proof.lookup_grand_product_commitment.point_mul(factor);
            r.point_add_assign(scaled);
        }
        PairingsBn254.Fr memory z_in_domain_size = state.z.pow(vk.domain_size);
        PairingsBn254.G1Point memory quotient_result = proof.quotient_poly_parts_commitments[0];
        PairingsBn254.Fr memory current_z = z_in_domain_size;
        PairingsBn254.G1Point memory tp;
        for(uint256 i = 1; i < proof.quotient_poly_parts_commitments.length; i++){
            tp = proof.quotient_poly_parts_commitments[i];
            tp.point_mul_assign(current_z);
            quotient_result.point_add_assign(tp);

            current_z.mul_assign(z_in_domain_size);
        }

        // TODO allocate constant sized array during code gen
        PairingsBn254.G1Point[2*SW-1 + 3] memory commitments_queried_at_z;
        PairingsBn254.Fr[2*SW-1 + 3] memory values_queried_at_z;
        uint256 idx = 0;
        commitments_queried_at_z[idx] = quotient_result;
        values_queried_at_z[idx] = proof.quotient_poly_opening_at_z;
        idx+=1;
        commitments_queried_at_z[idx] = r;
        values_queried_at_z[idx] = proof.linearization_poly_opening_at_z;
        for(uint256 i = 0; i<SW; i++){
            commitments_queried_at_z[idx] = proof.state_polys_commitments[i];
            values_queried_at_z[idx] = proof.state_polys_openings_at_z[i];
            idx+=1;
        }
        require(proof.gate_selectors_openings_at_z.length == 1);
        commitments_queried_at_z[idx] = vk.gate_selectors_commitments[0];
        values_queried_at_z[idx] = proof.gate_selectors_openings_at_z[0];
        idx+=1;
        for(uint256 i = 0; i<SW-1; i++){
            commitments_queried_at_z[idx] = vk.permutation_commitments[i];
            values_queried_at_z[idx] = proof.copy_permutation_polys_openings_at_z[i];
            idx+=1;
        }

        PairingsBn254.G1Point[2 + 3] memory commitments_queried_at_z_omega;
        PairingsBn254.Fr[2 + 3] memory values_queried_at_z_omega;

        commitments_queried_at_z_omega[0] = proof.copy_permutation_grand_product_commitment;
        commitments_queried_at_z_omega[1] = proof.state_polys_commitments[SW-1];

        values_queried_at_z_omega[0] = proof.copy_permutation_grand_product_opening_at_z_omega;
        values_queried_at_z_omega[1] = proof.state_polys_openings_at_z_omega[0];

        if(has_lookup){
            PairingsBn254.G1Point memory lookup_t_poly_commitment_aggregated = vk.lookup_tables_commitments[0];
            PairingsBn254.Fr memory current_eta = state.eta;
            for(uint256 i = 1; i < vk.lookup_tables_commitments.length; i++){
                tp = vk.lookup_tables_commitments[i].point_mul(current_eta);
                lookup_t_poly_commitment_aggregated.point_add_assign(tp);

                current_eta.mul_assign(state.eta);
            }
            commitments_queried_at_z[idx] = lookup_t_poly_commitment_aggregated;
            values_queried_at_z[idx] = proof.lookup_t_poly_opening_at_z;
            idx +=1;
            commitments_queried_at_z[idx] = vk.lookup_selector_commitment;
            values_queried_at_z[idx] = proof.lookup_selector_poly_opening_at_z;
            idx +=1;
            commitments_queried_at_z[idx] = vk.lookup_table_type_commitment;
            values_queried_at_z[idx] = proof.lookup_table_type_poly_opening_at_z;

            commitments_queried_at_z_omega[2] = proof.lookup_s_poly_commitment;
            values_queried_at_z_omega[2] = proof.lookup_s_poly_opening_at_z_omega;
            commitments_queried_at_z_omega[3] = proof.lookup_grand_product_commitment;
            values_queried_at_z_omega[3] = proof.lookup_grand_product_opening_at_z_omega;
            commitments_queried_at_z_omega[4] = lookup_t_poly_commitment_aggregated;
            values_queried_at_z_omega[4] = proof.lookup_t_poly_opening_at_z_omega;
        }

        require(commitments_queried_at_z_omega.length == values_queried_at_z_omega.length);

        idx = 0;
        PairingsBn254.G1Point memory  aggregated_commitment_at_z = commitments_queried_at_z[idx];
        PairingsBn254.Fr memory  aggregated_opening_at_z = values_queried_at_z[idx];
        idx+=1;
        PairingsBn254.Fr memory  aggregation_challenge = one;

        for(uint256 i = idx; i < commitments_queried_at_z.length; i++){
            aggregation_challenge.mul_assign(state.v);

            scaled = commitments_queried_at_z[idx].point_mul(aggregation_challenge);
            aggregated_commitment_at_z.point_add_assign(scaled);

            t = values_queried_at_z[idx];
            t.mul_assign(aggregation_challenge);
            aggregated_opening_at_z.add_assign(t);
            idx+=1;
        }

        aggregation_challenge.mul_assign(state.v);

        idx = 0;
        PairingsBn254.G1Point memory  aggregated_commitment_at_z_omega = commitments_queried_at_z_omega[idx];
        PairingsBn254.Fr memory  aggregated_opening_at_z_omega = values_queried_at_z_omega[idx];
        idx+=1;
        for(uint256 i = idx; i < commitments_queried_at_z_omega.length; i++){
            aggregation_challenge.mul_assign(state.v);

            scaled = commitments_queried_at_z_omega[idx].point_mul(aggregation_challenge);
            aggregated_commitment_at_z_omega.point_add_assign(scaled);

            t = values_queried_at_z_omega[idx];
            t.mul_assign(aggregation_challenge);
            aggregated_opening_at_z_omega.add_assign(t);
            idx+=1;
        }


        // q(x) = f(x) - f(z) / (x - z)
        // q(x) * (x-z)  = f(x) - f(z)
        // e(q(x), (x-z)) * e(f(x) - f(z), -G2) == 1

        // f(x)
        PairingsBn254.G1Point memory  pair_with_generator = aggregated_commitment_at_z;
        aggregated_commitment_at_z_omega.point_mul_assign(state.u);
        pair_with_generator.point_add_assign(aggregated_commitment_at_z_omega);

        // - f(z)*g
        PairingsBn254.Fr memory  aggregated_value = aggregated_opening_at_z_omega;
        aggregated_value.mul_assign(state.u);
        aggregated_value.add_assign(aggregated_opening_at_z);
        tp = PairingsBn254.P1().point_mul(aggregated_value);
        pair_with_generator.point_sub_assign(tp);

        // +z * q(x)
        tp = proof.opening_proof_at_z.point_mul(state.z);
        t = z_omega;
        t.mul_assign(state.u);
        PairingsBn254.G1Point memory t1 = proof.opening_proof_at_z_omega.point_mul(t);
        tp.point_add_assign(t1);

        // rhs
        PairingsBn254.G1Point memory pair_with_x = proof.opening_proof_at_z_omega.point_mul(state.u);
        pair_with_x.point_add_assign(proof.opening_proof_at_z);
        pair_with_x.negate();

        bool valid = PairingsBn254.pairingProd2(pair_with_generator, vk.g2_elements[0], pair_with_x, vk.g2_elements[1]);
        
        return valid;
    }
}