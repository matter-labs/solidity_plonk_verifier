pragma solidity >=0.5.0 <0.7.0;

library PairingsBn254 {
    uint256 constant q_mod = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
    uint256 constant r_mod = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
    
    struct G1Point {
        uint256 X;
        uint256 Y;
    } 
    
    struct Fr {
        uint256 value;
    }
    
    function newFr(uint256 fr) internal pure returns (Fr memory) {
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
        assert(fr.value != 0);
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
    
    function newG1(uint256 x, uint256 y) internal pure returns (G1Point memory) {
        return G1Point(x, y);
    }
    
    function newG2(uint256[2] memory x, uint256[2] memory y) internal pure returns (G2Point memory) {
        return G2Point(x, y);
    }
    
    function copyG1(G1Point memory self) internal pure returns (G1Point memory result) {
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
        if (self.X == 0 && self.Y == 0)
            return;
        self.Y = q_mod - self.Y;
    }
    
    function is_infinity(G1Point memory p) internal pure returns (bool) {
        if (p.X == 0 && p.Y == 0) {
            return true;
        }
        return false;
    }

    function pointAdd(G1Point memory p1, G1Point memory p2)
        internal view returns (G1Point memory r)
    {
        uint256[4] memory input;
        input[0] = p1.X;
        input[1] = p1.Y;
        input[2] = p2.X;
        input[3] = p2.Y;
        bool success = false;
        assembly {
            success := staticcall(gas(), 6, input, 0x80, r, 0x40)
        }
        require(success);
    }
    
    function point_add_assign(G1Point memory p1, G1Point memory p2)
        internal view
    {
        uint256[4] memory input;
        input[0] = p1.X;
        input[1] = p1.Y;
        input[2] = p2.X;
        input[3] = p2.Y;
        bool success = false;
        assembly {
            success := staticcall(gas(), 6, input, 0x80, p1, 0x40)
        }
        require(success);
    }
    
    function point_sub_assign(G1Point memory p1, G1Point memory p2)
        internal view
    {
        uint256[4] memory input;
        if (p1.X == 0 && p1.Y == 0) {
            p1.X = p2.X;
            p1.Y = q_mod - p2.Y;
            return;
        } else if (p2.X == 0 && p2.Y == 0) {
            return;
        } else {
            input[0] = p1.X;
            input[1] = p1.Y;
            input[2] = p2.X;
            input[3] = q_mod - p2.Y;
        }
        bool success = false;
        assembly {
            success := staticcall(gas(), 6, input, 0x80, p1, 0x40)
        }
        require(success);
    }


    function pointMul(G1Point memory p, Fr memory s)
        internal view returns (G1Point memory r)
    {
        uint[3] memory input;
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s.value;
        bool success;
        assembly {
            success := staticcall(gas(), 7, input, 0x60, r, 0x40)
        }
        require(success);
    }
    
    function point_mul_assign(G1Point memory p, Fr memory s)
        internal view
    {
        uint[3] memory input;
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s.value;
        bool success;
        assembly {
            success := staticcall(gas(), 7, input, 0x60, p, 0x40)
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
        self.state_0 = keccak256(abi.encodePacked(DST_0, self.state_0, value));
        self.state_1 = keccak256(abi.encodePacked(DST_1, self.state_1, value));
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

contract Plonk4VerifierWithAccessToDNext {
    using PairingsBn254 for PairingsBn254.G1Point;
    using PairingsBn254 for PairingsBn254.G2Point;
    using PairingsBn254 for PairingsBn254.Fr;
    
    using TranscriptLibrary for TranscriptLibrary.Transcript;
    
    struct VerificationKey {
        uint256 domain_size;
        uint256 num_inputs;
        PairingsBn254.Fr omega;
        PairingsBn254.G1Point[] selector_commitments;
        PairingsBn254.G1Point[] next_step_selector_commitments;
        PairingsBn254.G1Point[] permutation_commitments;
        PairingsBn254.Fr[] permutation_non_residues;
        PairingsBn254.G2Point g2_x;
    }
    
    struct Proof {
        uint256[] input_values;
        PairingsBn254.G1Point[] wire_commitments;
        PairingsBn254.G1Point grand_product_commitment;
        PairingsBn254.G1Point[] quotient_poly_commitments;
        PairingsBn254.Fr[] wire_values_at_z;
        PairingsBn254.Fr[] wire_values_at_z_omega;
        PairingsBn254.Fr grand_product_at_z_omega;
        PairingsBn254.Fr quotient_polynomial_at_z;
        PairingsBn254.Fr linearization_polynomial_at_z;
        PairingsBn254.Fr[] permutation_polynomials_at_z;
    
        PairingsBn254.G1Point opening_at_z_proof;
        PairingsBn254.G1Point opening_at_z_omega_proof;
    }
    
    struct PartialVerifierState {
        PairingsBn254.Fr alpha;
        PairingsBn254.Fr beta;
        PairingsBn254.Fr gamma;
        PairingsBn254.Fr v;
        PairingsBn254.Fr u;
        PairingsBn254.Fr z;
        PairingsBn254.Fr[] cached_lagrange_evals;
    }
    
    uint256 constant STATE_WIDTH = 4;
    
    function evaluate_lagrange_poly_out_of_domain(
        uint256 poly_num, 
        uint256 domain_size, 
        PairingsBn254.Fr memory omega, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr memory res) {
        require(poly_num < domain_size);
        PairingsBn254.Fr memory one = PairingsBn254.newFr(1);
        PairingsBn254.Fr memory omega_power = omega.pow(poly_num);
        res = at.pow(domain_size);
        res.sub_assign(one);
        assert(res.value != 0); // Vanishing polynomial can not be zero at point `at`
        res.mul_assign(omega_power);
        
        PairingsBn254.Fr memory den = PairingsBn254.copy(at);
        den.sub_assign(omega_power);
        den.mul_assign(PairingsBn254.newFr(domain_size));
        
        den = den.inverse();
        
        res.mul_assign(den);
    }
    
    function batch_evaluate_lagrange_poly_out_of_domain(
        uint256[] memory poly_nums, 
        uint256 domain_size, 
        PairingsBn254.Fr memory omega, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr[] memory res) {
        PairingsBn254.Fr memory one = PairingsBn254.newFr(1);
        PairingsBn254.Fr memory tmp_1 = PairingsBn254.newFr(0);
        PairingsBn254.Fr memory tmp_2 = PairingsBn254.newFr(domain_size);
        PairingsBn254.Fr memory vanishing_at_z = at.pow(domain_size);
        vanishing_at_z.sub_assign(one);
        // we can not have random point z be in domain
        assert(vanishing_at_z.value != 0);
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
        partial_products[0].assign(PairingsBn254.newFr(1));
        for (uint i = 1; i < dens.length - 1; i++) {
            partial_products[i].assign(dens[i-1]);
            partial_products[i].mul_assign(dens[i]);
        }
    
        tmp_2.assign(partial_products[partial_products.length - 1]);
        tmp_2.mul_assign(dens[dens.length - 1]);
        tmp_2 = tmp_2.inverse(); // tmp_2 contains a^-1 * b^-1 (with! the last one)
        
        for (uint i = dens.length - 1; i < dens.length; i--) {
            dens[i].assign(tmp_2); // all inversed
            dens[i].mul_assign(partial_products[i]); // clear lowest terms    
            tmp_2.mul_assign(dens[i]);
        }
        
        for (uint i = 0; i < nums.length; i--) {
            nums[i].mul_assign(dens[i]);
        }

        return nums;
    }
    
    function evaluate_vanishing(
        uint256 domain_size, 
        PairingsBn254.Fr memory at
    ) internal view returns (PairingsBn254.Fr memory res) {
        res = at.pow(domain_size);
        res.sub_assign(PairingsBn254.newFr(1));
    }
    
    function verify_at_z(
        PartialVerifierState memory state,
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (bool) {
        PairingsBn254.Fr memory lhs = evaluate_vanishing(vk.domain_size, state.z);
        assert(lhs.value != 0); // we can not check a polynomial relationship if point `z` is in the domain
        lhs.mul_assign(proof.quotient_polynomial_at_z);
    
        PairingsBn254.Fr memory quotient_challenge = PairingsBn254.newFr(1);
        PairingsBn254.Fr memory rhs = PairingsBn254.copy(proof.linearization_polynomial_at_z);
        
        // public inputs
        PairingsBn254.Fr memory tmp = PairingsBn254.newFr(0);
        for (uint256 i = 0; i < proof.input_values.length; i++) {
            tmp.assign(state.cached_lagrange_evals[i]);
            tmp.mul_assign(PairingsBn254.newFr(proof.input_values[i]));
            rhs.add_assign(tmp);
        }
        
        quotient_challenge.mul_assign(state.alpha);
        
        PairingsBn254.Fr memory z_part = PairingsBn254.copy(proof.grand_product_at_z_omega);
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            tmp.assign(proof.permutation_polynomials_at_z[i]);
            tmp.mul_assign(state.beta);
            tmp.add_assign(state.gamma);
            tmp.add_assign(proof.wire_values_at_z[i]);
            
            z_part.mul_assign(tmp);
        }
        
        tmp.assign(state.gamma);
        // we need a wire value of the last polynomial in enumeration
        tmp.add_assign(proof.wire_values_at_z[STATE_WIDTH - 1]);
        
        z_part.mul_assign(tmp);
        z_part.mul_assign(quotient_challenge);
        
        rhs.sub_assign(z_part);
        
        quotient_challenge.mul_assign(state.alpha);
        
        tmp.assign(state.cached_lagrange_evals[0]);
        tmp.mul_assign(quotient_challenge);
        
        rhs.sub_assign(tmp);
        
        return lhs.value == rhs.value;
    }
    
    function reconstruct_d(
        PartialVerifierState memory state,
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (PairingsBn254.G1Point memory res) {
        // we compute what power of v is used as a delinearization factor in batch opening of 
        // commitments. Let's label W(x) = 1 / (x - z) *
        // [
            // t_0(x) + z^n * t_1(x) + z^2n * t_2(x) + z^3n * t_3(x) - t(z)
            // + v (r(x) - r(z))
            // + v^{2..5} * (witness(x) - witness(z))
            // + v^(6..8) * (permutation(x) - permutation(z))
        // ]
        // W'(x) = 1 / (x - z*omega) *
        // [
            // + v^9 (z(x) - z(z*omega)) <- we need this power
            // + v^10 * (d(x) - d(z*omega))
        // ]
        //
        // we pay a little for a few arithmetic operations to not introduce another constant
        uint256 power_for_z_omega_opening = 1 + 1 + STATE_WIDTH + STATE_WIDTH - 1;
        res = PairingsBn254.copyG1(vk.selector_commitments[STATE_WIDTH + 1]);
                
        PairingsBn254.G1Point memory tmp_g1 = PairingsBn254.P1();
        PairingsBn254.Fr memory tmp_fr = PairingsBn254.newFr(0);
        
        // addition gates
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            tmp_g1 = vk.selector_commitments[i].pointMul(proof.wire_values_at_z[i]);
            res.point_add_assign(tmp_g1);
        }
        
        // multiplication gate
        tmp_fr.assign(proof.wire_values_at_z[0]);
        tmp_fr.mul_assign(proof.wire_values_at_z[1]);
        tmp_g1 = vk.selector_commitments[STATE_WIDTH].pointMul(tmp_fr);
        res.point_add_assign(tmp_g1);
        
        // d_next
        tmp_g1 = vk.next_step_selector_commitments[0].pointMul(proof.wire_values_at_z_omega[0]);
        res.point_add_assign(tmp_g1);
        
        // z * non_res * beta + gamma + a
        PairingsBn254.Fr memory grand_product_part_at_z = PairingsBn254.copy(state.z);
        grand_product_part_at_z.mul_assign(state.beta);
        grand_product_part_at_z.add_assign(proof.wire_values_at_z[0]);
        grand_product_part_at_z.add_assign(state.gamma);
        for (uint256 i = 0; i < vk.permutation_non_residues.length; i++) {
            tmp_fr.assign(state.z);
            tmp_fr.mul_assign(vk.permutation_non_residues[i]);
            tmp_fr.mul_assign(state.beta);
            tmp_fr.add_assign(state.gamma);
            tmp_fr.add_assign(proof.wire_values_at_z[i+1]);
            
            grand_product_part_at_z.mul_assign(tmp_fr);
        }
        
        grand_product_part_at_z.mul_assign(state.alpha);
    
        tmp_fr.assign(state.cached_lagrange_evals[0]);
        tmp_fr.mul_assign(state.alpha);
        tmp_fr.mul_assign(state.alpha);
        
        grand_product_part_at_z.add_assign(tmp_fr);
        
        PairingsBn254.Fr memory grand_product_part_at_z_omega = state.v.pow(power_for_z_omega_opening);
        grand_product_part_at_z_omega.mul_assign(state.u);
        
        PairingsBn254.Fr memory last_permutation_part_at_z = PairingsBn254.newFr(1);
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            tmp_fr.assign(state.beta);
            tmp_fr.mul_assign(proof.permutation_polynomials_at_z[i]);
            tmp_fr.add_assign(state.gamma);
            tmp_fr.add_assign(proof.wire_values_at_z[i]);
            
            last_permutation_part_at_z.mul_assign(tmp_fr);
        }

        last_permutation_part_at_z.mul_assign(state.beta);
        last_permutation_part_at_z.mul_assign(proof.grand_product_at_z_omega);
        last_permutation_part_at_z.mul_assign(state.alpha);
        
        // add to the linearization
        tmp_g1 = proof.grand_product_commitment.pointMul(grand_product_part_at_z);
        tmp_g1.point_sub_assign(vk.permutation_commitments[STATE_WIDTH - 1].pointMul(last_permutation_part_at_z));

        res.point_add_assign(tmp_g1);
        res.point_mul_assign(state.v);
        
        res.point_add_assign(proof.grand_product_commitment.pointMul(grand_product_part_at_z_omega));
    }
    
    function verify_commitments(
        PartialVerifierState memory state,
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (bool) {
        PairingsBn254.G1Point memory d = reconstruct_d(state, proof, vk);
        
        PairingsBn254.Fr memory z_in_domain_size = state.z.pow(vk.domain_size);
        
        PairingsBn254.G1Point memory tmp_g1 = PairingsBn254.P1();

        PairingsBn254.Fr memory aggregation_challenge = PairingsBn254.newFr(1);
        
        PairingsBn254.G1Point memory commitment_aggregation = PairingsBn254.copyG1(proof.quotient_poly_commitments[0]);
        PairingsBn254.Fr memory tmp_fr = PairingsBn254.newFr(1);
        for (uint i = 1; i < proof.quotient_poly_commitments.length; i++) {
            tmp_fr.mul_assign(z_in_domain_size);
            tmp_g1 = proof.quotient_poly_commitments[i].pointMul(tmp_fr);
            commitment_aggregation.point_add_assign(tmp_g1);
        }

        aggregation_challenge.mul_assign(state.v);
        commitment_aggregation.point_add_assign(d);
        
        for (uint i = 0; i < proof.wire_commitments.length; i++) {
            aggregation_challenge.mul_assign(state.v);
            tmp_g1 = proof.wire_commitments[i].pointMul(aggregation_challenge);
            commitment_aggregation.point_add_assign(tmp_g1);
        }
        
        for (uint i = 0; i < vk.permutation_commitments.length - 1; i++) {
            aggregation_challenge.mul_assign(state.v);
            tmp_g1 = vk.permutation_commitments[i].pointMul(aggregation_challenge);
            commitment_aggregation.point_add_assign(tmp_g1);
        }
        
        aggregation_challenge.mul_assign(state.v);

        aggregation_challenge.mul_assign(state.v);

        tmp_fr.assign(aggregation_challenge);
        tmp_fr.mul_assign(state.u);
        tmp_g1 = proof.wire_commitments[STATE_WIDTH - 1].pointMul(tmp_fr);
        commitment_aggregation.point_add_assign(tmp_g1);
        
        // collect opening values
        aggregation_challenge = PairingsBn254.newFr(1);
        
        PairingsBn254.Fr memory aggregated_value = PairingsBn254.copy(proof.quotient_polynomial_at_z);
        
        aggregation_challenge.mul_assign(state.v);

        tmp_fr.assign(proof.linearization_polynomial_at_z);
        tmp_fr.mul_assign(aggregation_challenge);
        aggregated_value.add_assign(tmp_fr);
        
        for (uint i = 0; i < proof.wire_values_at_z.length; i++) {
            aggregation_challenge.mul_assign(state.v);
            
            tmp_fr.assign(proof.wire_values_at_z[i]);
            tmp_fr.mul_assign(aggregation_challenge);
            aggregated_value.add_assign(tmp_fr);
        }
        
        for (uint i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            aggregation_challenge.mul_assign(state.v);

            tmp_fr.assign(proof.permutation_polynomials_at_z[i]);
            tmp_fr.mul_assign(aggregation_challenge);
            aggregated_value.add_assign(tmp_fr);
        }
        
        aggregation_challenge.mul_assign(state.v);

        tmp_fr.assign(proof.grand_product_at_z_omega);
        tmp_fr.mul_assign(aggregation_challenge);
        tmp_fr.mul_assign(state.u);
        aggregated_value.add_assign(tmp_fr);
        
        aggregation_challenge.mul_assign(state.v);

        tmp_fr.assign(proof.wire_values_at_z_omega[0]);
        tmp_fr.mul_assign(aggregation_challenge);
        tmp_fr.mul_assign(state.u);
        aggregated_value.add_assign(tmp_fr);
        
        commitment_aggregation.point_sub_assign(PairingsBn254.P1().pointMul(aggregated_value));
        
        PairingsBn254.G1Point memory pair_with_generator = commitment_aggregation;
        pair_with_generator.point_add_assign(proof.opening_at_z_proof.pointMul(state.z));
        
        tmp_fr.assign(state.z);
        tmp_fr.mul_assign(vk.omega);
        tmp_fr.mul_assign(state.u);
        pair_with_generator.point_add_assign(proof.opening_at_z_omega_proof.pointMul(tmp_fr));
        
        PairingsBn254.G1Point memory pair_with_x = proof.opening_at_z_omega_proof.pointMul(state.u);
        pair_with_x.point_add_assign(proof.opening_at_z_proof);
        pair_with_x.negate();
        
        return PairingsBn254.pairingProd2(pair_with_generator, PairingsBn254.P2(), pair_with_x, vk.g2_x);
    }
    
    function verify_initial(
        PartialVerifierState memory state, 
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (bool) {
        require(proof.input_values.length == vk.num_inputs);
        require(vk.num_inputs >= 1);
        TranscriptLibrary.Transcript memory transcript = TranscriptLibrary.new_transcript();
        for (uint256 i = 0; i < vk.num_inputs; i++) {
            transcript.update_with_u256(proof.input_values[i]);
        }
        
        for (uint256 i = 0; i < proof.wire_commitments.length; i++) {
            transcript.update_with_g1(proof.wire_commitments[i]);
        }
        
        state.beta = transcript.get_challenge();
        state.gamma = transcript.get_challenge();
        
        transcript.update_with_g1(proof.grand_product_commitment);
        state.alpha = transcript.get_challenge();
        
        for (uint256 i = 0; i < proof.quotient_poly_commitments.length; i++) {
            transcript.update_with_g1(proof.quotient_poly_commitments[i]);
        }
    
        state.z = transcript.get_challenge();
        
        uint256[] memory lagrange_poly_numbers = new uint256[](vk.num_inputs);
        for (uint256 i = 0; i < lagrange_poly_numbers.length; i++) {
            lagrange_poly_numbers[i] = i;
        }
        
        state.cached_lagrange_evals = batch_evaluate_lagrange_poly_out_of_domain(
            lagrange_poly_numbers,
            vk.domain_size, 
            vk.omega, state.z
        );

        bool valid = verify_at_z(state, proof, vk);

        if (valid == false) {
            return false;
        }
        
        for (uint256 i = 0; i < proof.wire_values_at_z.length; i++) {
            transcript.update_with_fr(proof.wire_values_at_z[i]);
        }
        
        for (uint256 i = 0; i < proof.wire_values_at_z_omega.length; i++) {
            transcript.update_with_fr(proof.wire_values_at_z_omega[i]);
        }
        
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            transcript.update_with_fr(proof.permutation_polynomials_at_z[i]);
        }
        
        transcript.update_with_fr(proof.quotient_polynomial_at_z);
        transcript.update_with_fr(proof.linearization_polynomial_at_z);
        
        state.v = transcript.get_challenge();
        transcript.update_with_g1(proof.opening_at_z_proof);
        transcript.update_with_g1(proof.opening_at_z_omega_proof);
        state.u = transcript.get_challenge();
        
        return true;
    }

    // This verifier is for a PLONK with a state width 4
    // and main gate equation
    // q_a(X) * a(X) + 
    // q_b(X) * b(X) + 
    // q_c(X) * c(X) + 
    // q_d(X) * d(X) + 
    // q_m(X) * a(X) * b(X) + 
    // q_constants(X)+ 
    // q_d_next(X) * d(X*omega)
    // where q_{}(X) are selectors a, b, c, d - state (witness) polynomials
    // q_d_next(X) "peeks" into the next row of the trace, so it takes 
    // the same d(X) polynomial, but shifted  
    
    function verify(Proof memory proof, VerificationKey memory vk) internal view returns (bool) {
        PartialVerifierState memory state;
        
        bool valid = verify_initial(state, proof, vk);
        
        if (valid == false) {
            return false;
        }
        
        valid = verify_commitments(state, proof, vk);
        
        return valid;
    }
}

contract ConcreteVerifier is Plonk4VerifierWithAccessToDNext {
    function get_verification_key() internal pure returns(VerificationKey memory vk) {
        vk.domain_size = {{domain_size}};
        vk.num_inputs = {{num_inputs}};
        vk.omega = PairingsBn254.newFr({{omega}});
        vk.selector_commitments = new PairingsBn254.G1Point[](STATE_WIDTH+2);
        vk.selector_commitments[0] = PairingsBn254.newG1(
            {{selector_commitment_0_0}},
            {{selector_commitment_0_1}}
        );
        vk.selector_commitments[1] = PairingsBn254.newG1(
            {{selector_commitment_1_0}},
            {{selector_commitment_1_1}}
        );
        vk.selector_commitments[2] = PairingsBn254.newG1(
            {{selector_commitment_2_0}},
            {{selector_commitment_2_1}}
        );
        vk.selector_commitments[3] = PairingsBn254.newG1(
            {{selector_commitment_3_0}},
            {{selector_commitment_3_1}}
        );
        vk.selector_commitments[4] = PairingsBn254.newG1(
            {{selector_commitment_4_0}},
            {{selector_commitment_4_1}}
        );
        vk.selector_commitments[5] = PairingsBn254.newG1(
            {{selector_commitment_5_0}},
            {{selector_commitment_5_1}}
        );
        
        // we only have access to value of the d(x) witness polynomial on the next
        // trace step, so we only need one element here and deal with it in other places
        // by having this in mind
        vk.next_step_selector_commitments = new PairingsBn254.G1Point[](1);
        vk.next_step_selector_commitments[0] = PairingsBn254.newG1(
            {{next_step_selector_commitment_0_0}},
            {{next_step_selector_commitment_0_1}}
        );
        
        vk.permutation_commitments = new PairingsBn254.G1Point[](STATE_WIDTH);
         vk.permutation_commitments[0] = PairingsBn254.newG1(
            {{permutation_commitment_0_0}},
            {{permutation_commitment_0_1}}
        );
        vk.permutation_commitments[1] = PairingsBn254.newG1(
            {{permutation_commitment_1_0}},
            {{permutation_commitment_1_1}}
        );
        vk.permutation_commitments[2] = PairingsBn254.newG1(
            {{permutation_commitment_2_0}},
            {{permutation_commitment_2_1}}
        );
        vk.permutation_commitments[3] = PairingsBn254.newG1(
            {{permutation_commitment_3_0}},
            {{permutation_commitment_3_1}}
        );
        
        vk.permutation_non_residues = new PairingsBn254.Fr[](STATE_WIDTH-1);
        vk.permutation_non_residues[0] = PairingsBn254.newFr(
            {{permutation_non_residue_0}}
        );
        vk.permutation_non_residues[1] = PairingsBn254.newFr(
            {{permutation_non_residue_1}}
        );
        vk.permutation_non_residues[2] = PairingsBn254.newFr(
            {{permutation_non_residue_2}}
        );
        
        vk.g2_x = PairingsBn254.newG2(
            [{{g2_x_x_c1}},
             {{g2_x_x_c0}}],
            [{{g2_x_y_c1}},
             {{g2_x_y_c0}}]
        );
    }


    function deserialize_proof(
        uint256 expected_inputs, 
        uint256[] memory public_inputs, 
        uint256[] memory serialized_proof
    ) internal pure returns(Proof memory proof) {
        assert(expected_inputs == public_inputs.length);
        proof.input_values = new uint256[](expected_inputs);
        for (uint256 i = 0; i < expected_inputs; i++) {
            proof.input_values[i] = public_inputs[i];
        }
 
        uint256 j = 0;
        proof.wire_commitments = new PairingsBn254.G1Point[](STATE_WIDTH);
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.wire_commitments[i] = PairingsBn254.newG1(
                serialized_proof[j],
                serialized_proof[j+1]
            );

            j += 2;
        }
        
        proof.grand_product_commitment = PairingsBn254.newG1(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;
        
        proof.quotient_poly_commitments = new PairingsBn254.G1Point[](STATE_WIDTH);
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.quotient_poly_commitments[i] = PairingsBn254.newG1(
                serialized_proof[j],
                serialized_proof[j+1]
            );

            j += 2;
        }
        
        proof.wire_values_at_z = new PairingsBn254.Fr[](STATE_WIDTH);
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.wire_values_at_z[i] = PairingsBn254.newFr(
                serialized_proof[j]
            );

            j += 1;
        }
        
        proof.wire_values_at_z_omega = new PairingsBn254.Fr[](1);
        for (uint256 i = 0; i < proof.wire_values_at_z_omega.length; i++) {
            proof.wire_values_at_z_omega[i] = PairingsBn254.newFr(
                serialized_proof[j]
            );

            j += 1;
        }
        
        proof.grand_product_at_z_omega = PairingsBn254.newFr(
                serialized_proof[j]
            );

        j += 1;

        proof.quotient_polynomial_at_z = PairingsBn254.newFr(
            serialized_proof[j]
        );

        j += 1;

        proof.linearization_polynomial_at_z = PairingsBn254.newFr(
            serialized_proof[j]
        );

        j += 1;
    
        proof.permutation_polynomials_at_z = new PairingsBn254.Fr[](STATE_WIDTH - 1);
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            proof.permutation_polynomials_at_z[i] = PairingsBn254.newFr(
                serialized_proof[j]
            );

            j += 1;
        }

        proof.opening_at_z_proof = PairingsBn254.newG1(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;

        proof.opening_at_z_omega_proof = PairingsBn254.newG1(
                serialized_proof[j],
                serialized_proof[j+1]
        );
    
        j += 2;
        assert(j == serialized_proof.length);
    }
    
    function verify(
        uint256[] memory public_inputs, 
        uint256[] memory serialized_proof
    ) public view returns (bool) {
        VerificationKey memory vk = get_verification_key();
        uint256 expected_inputs = vk.num_inputs;

        Proof memory proof = deserialize_proof(expected_inputs, public_inputs, serialized_proof);

        bool valid = verify(proof, vk);

        return valid;
    }  
}

// contract TestVerifier is Plonk4VerifierWithAccessToDNext {
//     function get_verification_key() internal pure returns(VerificationKey memory vk) {
//         vk.domain_size = 8;
//         vk.num_inputs = 1;
//         vk.omega = PairingsBn254.newFr(0x2b337de1c8c14f22ec9b9e2f96afef3652627366f8170a0a948dad4ac1bd5e80);
//         vk.selector_commitments = new PairingsBn254.G1Point[](STATE_WIDTH+2);
//         vk.selector_commitments[0] = PairingsBn254.newG1(
//             0x221f0a4724496eef2d6f4a3feb4ffc0208de7bfdb737ea8e2b835e0616d061d2,
//             0x11c939ddcd22087b6fd59d799e15924573b8a815a789d237414b8430fcf69974
//         );
//         vk.selector_commitments[1] = PairingsBn254.newG1(
//             0x1a491b60cdae811dc03c00655a64073cf536115324b340895c6db7e387eea681,
//             0x06e63b071abae5ad0dc281f103e234f34535d795699dcc6a895d5964422a5952
//         );
//         vk.selector_commitments[2] = PairingsBn254.newG1(
//             0x1a491b60cdae811dc03c00655a64073cf536115324b340895c6db7e387eea681,
//             0x297e136bc676ba7caa8dc3c57d9f236a524b92fbfed3fe22b2c332b29652a3f5
//         );
//         vk.selector_commitments[3] = PairingsBn254.newG1(
//             0x1a491b60cdae811dc03c00655a64073cf536115324b340895c6db7e387eea681,
//             0x297e136bc676ba7caa8dc3c57d9f236a524b92fbfed3fe22b2c332b29652a3f5
//         );
//         vk.selector_commitments[4] = PairingsBn254.newG1(
//             0x1e52d5ac827d81147c8e884f7b220e2779c7a8cd67df5ff5857995ebaa622597,
//             0x2b9fc3b60b69dee94225aef2b92a5bc716b180d76ab61cb8f37b04336f281022
//         );
//         vk.selector_commitments[5] = PairingsBn254.newG1(
//             0x0000000000000000000000000000000000000000000000000000000000000000,
//             0x0000000000000000000000000000000000000000000000000000000000000000
//         );
        
//         vk.next_step_selector_commitments = new PairingsBn254.G1Point[](1);
//         vk.next_step_selector_commitments[0] = PairingsBn254.newG1(
//             0x0000000000000000000000000000000000000000000000000000000000000000,
//             0x0000000000000000000000000000000000000000000000000000000000000000
//         );
        
//         vk.permutation_commitments = new PairingsBn254.G1Point[](STATE_WIDTH);
//          vk.permutation_commitments[0] = PairingsBn254.newG1(
//             0x055fac5e3320cc7dad48dd687d545df9ab5eaf2a5790798acf87aef48579f7ad,
//             0x1f555da7f16224290045aa71e459f34cb9eabbc4155fb68d6ebdab0af235dfdc
//         );
//         vk.permutation_commitments[1] = PairingsBn254.newG1(
//             0x0177ec222fce18b70e2e6f9c649612d73e6af69f24e8151b9cab607f37316084,
//             0x09d229cd3e9673a712806899bd09d65e7906dfd92b147bb656ab4bd713b8a8c8
//         );
//         vk.permutation_commitments[2] = PairingsBn254.newG1(
//             0x27c39bc902603bddf02576de7bfe7211889467877fa7b5611c495a0e462dc012,
//             0x093501f69f0617158db4f72c4f1183e284d2c1d92099add2713d48fc851bb698
//         );
//         vk.permutation_commitments[3] = PairingsBn254.newG1(
//             0x19c5d1b7df1125f24646c8e7ac90689c6cd593703e7cd4d7918e3d7c9c2f220a,
//             0x0daf49f7dfa3f8d741df0efcef1f171c961b47f8b28d46a2af00be85c8b6e713
//         );
        
//         vk.permutation_non_residues = new PairingsBn254.Fr[](STATE_WIDTH-1);
//         vk.permutation_non_residues[0] = PairingsBn254.newFr(
//             0x0000000000000000000000000000000000000000000000000000000000000005
//         );
//         vk.permutation_non_residues[1] = PairingsBn254.newFr(
//             0x0000000000000000000000000000000000000000000000000000000000000005
//         );
//         vk.permutation_non_residues[2] = PairingsBn254.newFr(
//             0x0000000000000000000000000000000000000000000000000000000000000005
//         );
        
//         vk.g2_x = PairingsBn254.newG2(
//             [0x12740934ba9615b77b6a49b06fcce83ce90d67b1d0e2a530069e3a7306569a91,
//              0x116da8c89a0d090f3d8644ada33a5f1c8013ba7204aeca62d66d931b99afe6e7],
//             [0x25222d9816e5f86b4a7dedd00d04acc5c979c18bd22b834ea8c6d07c0ba441db,
//              0x076441042e77b6309644b56251f059cf14befc72ac8a6157d30924e58dc4c172]
//         );
//     }
    
//     function get_dummy_proof() internal pure returns(Proof memory proof) {
//         proof.input_values = new uint256[](1);
//         proof.input_values[0] = 1;
//         proof.wire_commitments = new PairingsBn254.G1Point[](STATE_WIDTH);
//         proof.wire_commitments[0] = PairingsBn254.newG1(
//             0x1b5e58c31beca50e963688ba8d6be08b1c90e964135450d882876c8a5e8bc613,
//             0x2d2321db916c233e56596bba357913cdad1e962e9264f9c74773af210a5ecd29
//         );
//         proof.wire_commitments[1] = PairingsBn254.newG1(
//             0x2c4951b832621f1b7315181f70905ba8cd4865753d98da5c28a474092edadbda,
//             0x27b14de154fd0af1394e5ca8f8c64032ee91cd9e046cbd7a24aeed18db9c3334
//         );
//         proof.wire_commitments[2] = PairingsBn254.newG1(
//             0x1a491b60cdae811dc03c00655a64073cf536115324b340895c6db7e387eea681,
//             0x06e63b071abae5ad0dc281f103e234f34535d795699dcc6a895d5964422a5952
//         );
//         proof.wire_commitments[3] = PairingsBn254.newG1(
//             0x0000000000000000000000000000000000000000000000000000000000000000,
//             0x0000000000000000000000000000000000000000000000000000000000000000
//         );
        
//         proof.grand_product_commitment = PairingsBn254.newG1(
//             0x161b6ef60ada7b486c9bfbc23350365422fc3eff69aa0716c1319b696f41dc0b,
//             0x1cbc0c1bd268e96ac40a377067ee2ca383983870ce555db61ec2ab25fa3996c4
//         );
        
//         proof.quotient_poly_commitments = new PairingsBn254.G1Point[](STATE_WIDTH);
//         proof.quotient_poly_commitments[0] = PairingsBn254.newG1(
//             0x1d8f5d10f66ea913f4abd01d178e6c0367f6ed7e6ad86f0637eec170043b794c,
//             0x16339d848618af3adfc06827fe67a3b7c65a989b442d6e4266056005052bb72c
//         );
//         proof.quotient_poly_commitments[1] = PairingsBn254.newG1(
//             0x2012ea21d65fe7bd0be88491eba50e779becced0c72030babce4ea652fd5e92a,
//             0x28c63d4ef6bbceb1602d5e4df8c7fe8e33a2360241c81fc320f16910bb05ed09
//         );
//         proof.quotient_poly_commitments[2] = PairingsBn254.newG1(
//             0x0041c38c5baae54d3e9cdd56a85e88ed976f3b082bcfef2067658886372067ba,
//             0x062c2b76fa36f2bf0bcdaa17a27e8c10053bec7870c599086a6338d01dbd3262
//         );
//         proof.quotient_poly_commitments[3] = PairingsBn254.newG1(
//             0x0000000000000000000000000000000000000000000000000000000000000000,
//             0x0000000000000000000000000000000000000000000000000000000000000000
//         );
        
//         proof.wire_values_at_z = new PairingsBn254.Fr[](STATE_WIDTH);
//         proof.wire_values_at_z[0] = PairingsBn254.newFr(
//             0x1b6902ca32bc26a34f6b846eeb93ac69fa4581617da3b318460eca73dd2d18e0
//         );
//         proof.wire_values_at_z[1] = PairingsBn254.newFr(
//             0x2c9adab6d77eb011f3e3e511f2f0df1e8fc695601c7775f8dfa1d881b244dc4c
//         );
//         proof.wire_values_at_z[2] = PairingsBn254.newFr(
//             0x2a139912dcbb852a87c257a100793c63da16ca0fd146177e4f47fbd066706f2d
//         );
//         proof.wire_values_at_z[3] = PairingsBn254.newFr(
//             0x0000000000000000000000000000000000000000000000000000000000000000
//         );
        
//         proof.wire_values_at_z_omega = new PairingsBn254.Fr[](1);
//         proof.wire_values_at_z_omega[0] = PairingsBn254.newFr(
//             0x0000000000000000000000000000000000000000000000000000000000000000
//         );
        
//         proof.grand_product_at_z_omega = PairingsBn254.newFr(
//             0x09334ae6db6affb8c1efb113552128645b728f9b64a4b28466243a8fdf31dcc3
//         );
        
//         proof.quotient_polynomial_at_z = PairingsBn254.newFr(
//             0x06dcd12571bf33be988cfb99610eb0c3e6fd135cee07799afb177a9a710ac7b3
//         );
        
//         proof.linearization_polynomial_at_z = PairingsBn254.newFr(
//             0x013c7a8c3d10661074f62207bacc0d98105182b29a9b6f077c7561cf4aea90f7
//         );
        
//         proof.permutation_polynomials_at_z = new PairingsBn254.Fr[](STATE_WIDTH - 1);
//         proof.permutation_polynomials_at_z[0] = PairingsBn254.newFr(
//             0x1d8cdc7c9e8e4c7b07d69f62746da251c5de4f0e677f4b157621e47a097c0f97
//         );
//         proof.permutation_polynomials_at_z[1] = PairingsBn254.newFr(
//             0x284d14ca9b828c96398366151333fac8ca4af048e7a5bf424d5c77a2731a8132
//         );
//         proof.permutation_polynomials_at_z[2] = PairingsBn254.newFr(
//             0x00f766126ed12460d91886cccd2e5f467f1856aec6a8f39c3b732dd4883548ba
//         );
        
//         proof.opening_at_z_proof = PairingsBn254.newG1(
//             0x0f0ec46b791690bfcc86e7e98c65b431e9ff623bed4ba78ef5ba1ff9c9e24468,
//             0x23b2037afd719f15f230aabc2631d8dca79170d013cd0182608e36fb4f18da0f
//         );
//         proof.opening_at_z_omega_proof = PairingsBn254.newG1(
//             0x14cdb2e64f4deef40be7a70381246a752698e8b798d3f00f4dee0e2e216acc1d,
//             0x199a18224eb1774d3646008365ec59c333b385262d5a32e9131cdd1493ddeb7f
//         );
//     }

//     function test() public view returns (bool) {
//         bool valid = verify(get_dummy_proof(), get_verification_key());
//         assert(valid);
//         return valid;
//     }
// }