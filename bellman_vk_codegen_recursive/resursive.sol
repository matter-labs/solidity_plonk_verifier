pragma solidity >=0.5.0 <0.7.0;

library PairingsBn254 {
    uint256 constant q_mod = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
    uint256 constant r_mod = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
    uint256 constant bn254_b_coeff = 3;

    struct G1Point {
        uint256 X;
        uint256 Y;
    }

    function new_fr(uint256 fr) internal pure returns (uint256) {
        require(fr < r_mod, "F"); /* !ERROR: Element is not in the field */
        return fr;
    }

    function inverse(uint256 fr) internal view returns (uint256) {
        require(fr != 0, "MP"); /* !ERROR: Can't invert zero */
        return pow(fr, r_mod-2);
    }

    function addm(uint256 self, uint256 other) internal pure returns (uint256){
        return addmod(self, other, r_mod);
    }

    function subm(uint256 self, uint256 other) internal pure returns (uint256) {
        return addmod(self, r_mod - other, r_mod);
    }

    function mulm(uint256 self, uint256 other) internal pure returns (uint256) {
        return mulmod(self, other, r_mod);
    }

    function pow(uint256 self, uint256 power) internal view returns (uint256) {
        uint256[6] memory input = [32, 32, 32, self, power, r_mod];
        uint256[1] memory result;
        bool success;
        assembly {
            success := staticcall(gas(), 0x05, input, 0xc0, result, 0x20)
        }
        require(success, "XO"); /* !ERROR: Error raising to power */
        return result[0];
    }

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

    function new_g1_checked(uint256 x, uint256 y) internal pure returns (G1Point memory) {
        if (x == 0 && y == 0) {
            // point of infinity is (0,0)
            return G1Point(x, y);
        }

        // check encoding
        require(x < q_mod, "T4"); /* !ERROR: 'x' is not in the field */
        require(y < q_mod, "VL"); /* !ERROR: 'y' is not in the field */
        // check on curve
        uint256 lhs = mulmod(y, y, q_mod); // y^2
        uint256 rhs = mulmod(x, x, q_mod); // x^2
        rhs = mulmod(rhs, x, q_mod); // x^3
        rhs = addmod(rhs, bn254_b_coeff, q_mod); // x^3 + b
        require(lhs == rhs, "YZ"); /* !ERROR: Point is not on the curve */

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
            require(self.X == 0, "3A"); /* !ERROR: Invalid infinity point encoding */
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
            require(success, "5"); /* !ERROR: Error adding curve points */
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
            require(success, "TL"); /* !ERROR: Error subtracting curve points */
        }
    }

    function point_mul(G1Point memory p, uint256 s)
        internal view returns (G1Point memory r)
    {
        point_mul_into_dest(p, s, r);
        return r;
    }

    function point_mul_assign(G1Point memory p, uint256 s)
       internal view
    {
        point_mul_into_dest(p, s, p);
    }

    function point_mul_into_dest(G1Point memory p, uint256 s, G1Point memory dest)
        internal view
    {
        uint[3] memory input;
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s;
        bool success;
        assembly {
            success := staticcall(gas(), 7, input, 0x60, dest, 0x40)
        }
        require(success, "TG"); /* !ERROR: Error multiplying curve point */
    }

    function pairing(G1Point[] memory p1, G2Point[] memory p2)
        internal view returns (bool)
    {
        require(p1.length == p2.length, "N"); /* !ERROR: Point arrays differ in length */
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

        require(success, "Q"); /* !ERROR: Error computing pairing */
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

    function update_with_fr(Transcript memory self, uint256 value) internal pure {
        update_with_u256(self, value);
    }

    function update_with_g1(Transcript memory self, PairingsBn254.G1Point memory p) internal pure {
        update_with_u256(self, p.X);
        update_with_u256(self, p.Y);
    }

    function get_challenge(Transcript memory self) internal pure returns(uint256 challenge) {
        bytes32 query = keccak256(abi.encodePacked(DST_CHALLENGE, self.state_0, self.state_1, self.challenge_counter));
        self.challenge_counter += 1;
        challenge = uint256(query) & FR_MASK;
        return challenge;
    }
}

contract Plonk4VerifierWithAccessToDNext {
    using PairingsBn254 for PairingsBn254.G1Point;
    using PairingsBn254 for PairingsBn254.G2Point;
    using PairingsBn254 for uint256;

    using TranscriptLibrary for TranscriptLibrary.Transcript;

    uint256 constant STATE_WIDTH = 4;
    uint256 constant NUM_DIFFERENT_GATES = 2;
    uint256 constant NUM_SETUP_POLYS_FOR_MAIN_GATE = 7;
    uint256 constant NUM_SETUP_POLYS_RANGE_CHECK_GATE = 0;
    uint256 constant ACCESSIBLE_STATE_POLYS_ON_NEXT_STEP = 1;
    uint256 constant NUM_GATE_SELECTORS_OPENED_EXPLICITLY = 1;

    uint256 constant RECURSIVE_CIRCUIT_INPUT_COMMITMENT_MASK = 0x00ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff;
    uint256 constant LIMB_WIDTH = 68;
    uint256 constant MAX_LIMB_VALUE = 1 << LIMB_WIDTH;

    struct VerificationKey {
        uint256 domain_size;
        uint256 num_inputs;
        uint256 omega;
        PairingsBn254.G1Point[NUM_SETUP_POLYS_FOR_MAIN_GATE + NUM_SETUP_POLYS_RANGE_CHECK_GATE] gate_setup_commitments;
        PairingsBn254.G1Point[NUM_DIFFERENT_GATES] gate_selector_commitments;
        PairingsBn254.G1Point[STATE_WIDTH] copy_permutation_commitments;
        uint256[STATE_WIDTH-1] copy_permutation_non_residues;
        PairingsBn254.G2Point g2_x;
    }

    struct Proof {
        uint256[] input_values;
        PairingsBn254.G1Point[STATE_WIDTH] wire_commitments;
        PairingsBn254.G1Point copy_permutation_grand_product_commitment;
        PairingsBn254.G1Point[STATE_WIDTH] quotient_poly_commitments;
        uint256[STATE_WIDTH] wire_values_at_z;
        uint256[ACCESSIBLE_STATE_POLYS_ON_NEXT_STEP] wire_values_at_z_omega;
        uint256[NUM_GATE_SELECTORS_OPENED_EXPLICITLY] gate_selector_values_at_z;
        uint256 copy_grand_product_at_z_omega;
        uint256 quotient_polynomial_at_z;
        uint256 linearization_polynomial_at_z;
        uint256[STATE_WIDTH-1] permutation_polynomials_at_z;
        PairingsBn254.G1Point opening_at_z_proof;
        PairingsBn254.G1Point opening_at_z_omega_proof;
    }

    struct PartialVerifierState {
        uint256 alpha;
        uint256 beta;
        uint256 gamma;
        uint256 v;
        uint256 u;
        uint256 z;
        uint256[] cached_lagrange_evals;
    }

    function batch_evaluate_lagrange_poly_out_of_domain(
        uint256[] memory poly_nums,
        uint256 domain_size,
        uint256 omega,
        uint256 at
    ) internal view returns (uint256[] memory res) {
        uint256 tmp_1 = 0;
        uint256 tmp_2 = PairingsBn254.new_fr(domain_size);
        uint256 vanishing_at_z = at.pow(domain_size);
        vanishing_at_z = vanishing_at_z.subm(1);
        // we can not have random point z be in domain
        require(vanishing_at_z != 0, "PM"); /* !ERROR: Point z can't be in domain */
        uint256[] memory nums = new uint256[](poly_nums.length);
        uint256[] memory dens = new uint256[](poly_nums.length);
        // numerators in a form omega^i * (z^n - 1)
        // denoms in a form (z - omega^i) * N
        for (uint i = 0; i < poly_nums.length; i++) {
            tmp_1 = omega.pow(poly_nums[i]); // power of omega
            nums[i] = vanishing_at_z;
            nums[i] = nums[i].mulm(tmp_1);

            dens[i] = at; // (X - omega^i) * N
            dens[i] = dens[i].subm(tmp_1);
            dens[i] = dens[i].mulm(tmp_2); // mul by domain size
        }

        uint256[] memory partial_products = new uint256[](poly_nums.length);
        partial_products[0] = 1;
        for (uint i = 1; i < dens.length - 1; i++) {
            partial_products[i] = dens[i-1];
            partial_products[i] = partial_products[i].mulm(dens[i]);
        }

        tmp_2 = partial_products[partial_products.length - 1];
        tmp_2 = tmp_2.mulm(dens[dens.length - 1]);
        tmp_2 = tmp_2.inverse(); // tmp_2 contains a^-1 * b^-1 (with! the last one)

        for (uint i = dens.length - 1; i < dens.length; i--) {
            dens[i] = tmp_2; // all inversed
            dens[i] = dens[i].mulm(partial_products[i]); // clear lowest terms
            tmp_2 = tmp_2.mulm(dens[i]);
        }

        for (uint i = 0; i < nums.length; i++) {
            nums[i] = nums[i].mulm(dens[i]);
        }

        return nums;
    }

    function evaluate_vanishing(
        uint256 domain_size,
        uint256 at
    ) internal view returns (uint256 res) {
        res = at.pow(domain_size);
        res = res.subm(1);
        return res;
    }

    function verify_at_z(
        PartialVerifierState memory state,
        Proof memory proof,
        VerificationKey memory vk
    ) internal view returns (bool) {
        uint256 lhs = evaluate_vanishing(vk.domain_size, state.z);
        require(lhs != 0, "MG"); /* !ERROR: Can't have point z in the domain */ // we can not check a polynomial relationship if point `z` is in the domain
        lhs = lhs.mulm(proof.quotient_polynomial_at_z);

        uint256 quotient_challenge = 1;
        uint256 rhs = proof.linearization_polynomial_at_z;

        // public inputs
        uint256 tmp = 0;
        uint256 inputs_term = 0;
        for (uint256 i = 0; i < proof.input_values.length; i++) {
            tmp = state.cached_lagrange_evals[i];
            tmp = tmp.mulm(PairingsBn254.new_fr(proof.input_values[i]));
            inputs_term = inputs_term.addm(tmp);
        }

        inputs_term = inputs_term.mulm(proof.gate_selector_values_at_z[0]);
        rhs = rhs.addm(inputs_term);

        // now we need 5th power
        quotient_challenge = quotient_challenge.mulm(state.alpha);
        quotient_challenge = quotient_challenge.mulm(state.alpha);
        quotient_challenge = quotient_challenge.mulm(state.alpha);
        quotient_challenge = quotient_challenge.mulm(state.alpha);
        quotient_challenge = quotient_challenge.mulm(state.alpha);

        uint256 z_part = proof.copy_grand_product_at_z_omega;
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            tmp = proof.permutation_polynomials_at_z[i];
            tmp = tmp.mulm(state.beta);
            tmp = tmp.addm(state.gamma);
            tmp = tmp.addm(proof.wire_values_at_z[i]);

            z_part = z_part.mulm(tmp);
        }

        tmp = state.gamma;
        // we need a wire value of the last polynomial in enumeration
        tmp = tmp.addm( proof.wire_values_at_z[STATE_WIDTH - 1]);

        z_part = z_part.mulm(tmp);
        z_part = z_part.mulm(quotient_challenge);

        rhs = rhs.subm(z_part);

        quotient_challenge = quotient_challenge.mulm(state.alpha);

        tmp = state.cached_lagrange_evals[0];
        tmp = tmp.mulm(quotient_challenge);

        rhs = rhs.subm(tmp);

        return lhs == rhs;
    }

    function add_contribution_from_range_constraint_gates(
        PartialVerifierState memory state,
        Proof memory proof,
        uint256 current_alpha
    ) internal pure returns (uint256, uint256) {
        // now add contribution from range constraint gate
        // we multiply selector commitment by all the factors (alpha*(c - 4d)(c - 4d - 1)(..-2)(..-3) + alpha^2 * (4b - c)()()() + {} + {})

        uint256 res = 0;
        uint256 t0 = 0;
        uint256 t1 = 0;
        uint256 t2 = 0;

        for (uint256 i = 0; i < 3; i++) {
            current_alpha = current_alpha.mulm(state.alpha);

            // high - 4*low

            // this is 4*low
            t0 = proof.wire_values_at_z[3 - i];
            t0 = t0.mulm(4);

            // high
            t1 = proof.wire_values_at_z[2 - i];
            t1 = t1.subm(t0);

            // t0 is now t1 - {0,1,2,3}

            // first unroll manually for -0;
            t2 = t1;

            // -1
            t0 = t1;
            t0 = t0.subm(1);
            t2 = t2.mulm(t0);

            // -2
            t0 = t1;
            t0 = t0.subm(2);
            t2 = t2.mulm(t0);

            // -3
            t0 = t1;
            t0 = t0.subm(3);
            t2 = t2.mulm(t0);

            t2 = t2.mulm(current_alpha);

            res = res.addm(t2);
        }

        // now also d_next - 4a

        current_alpha = current_alpha.mulm(state.alpha);

        // high - 4*low

        // this is 4*low
        t0 = proof.wire_values_at_z[0];
        t0 = t0.mulm(4);

        // high
        t1 = proof.wire_values_at_z_omega[0];
        t1 = t1.subm(t0);

        // t0 is now t1 - {0,1,2,3}

        // first unroll manually for -0;
        t2 = t1;

        // -1
        t0 = t1;
        t0 = t0.subm(1);
        t2 = t2.mulm(t0);

        // -2
        t0 = t1;
        t0 = t0.subm(2);
        t2 = t2.mulm(t0);

        // -3
        t0 = t1;
        t0 = t0.subm(3);
        t2 = t2.mulm(t0);

        t2 = t2.mulm(current_alpha);

        res = res.addm(t2);

        return (current_alpha,res);
    }

    function reconstruct_linearization_commitment(
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
            // + v^{6} * (selector(x) - selector(z))
            // + v^{7..9} * (permutation(x) - permutation(z))
        // ]
        // W'(x) = 1 / (x - z*omega) *
        // [
            // + v^10 (z(x) - z(z*omega)) <- we need this power
            // + v^11 * (d(x) - d(z*omega))
        // ]
        //

        // we reconstruct linearization polynomial virtual selector
        // for that purpose we first linearize over main gate (over all it's selectors)
        // and multiply them by value(!) of the corresponding main gate selector
        res = PairingsBn254.copy_g1(vk.gate_setup_commitments[STATE_WIDTH + 1]); // index of q_const(x)

        PairingsBn254.G1Point memory tmp_g1 = PairingsBn254.P1();
        uint256 tmp_fr = 0;

        // addition gates
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            tmp_g1 = vk.gate_setup_commitments[i].point_mul(proof.wire_values_at_z[i]);
            res.point_add_assign(tmp_g1);
        }

        // multiplication gate
        tmp_fr = proof.wire_values_at_z[0];
        tmp_fr = tmp_fr.mulm(proof.wire_values_at_z[1]);
        tmp_g1 = vk.gate_setup_commitments[STATE_WIDTH].point_mul(tmp_fr);
        res.point_add_assign(tmp_g1);

        // d_next
        tmp_g1 = vk.gate_setup_commitments[STATE_WIDTH+2].point_mul(proof.wire_values_at_z_omega[0]); // index of q_d_next(x)
        res.point_add_assign(tmp_g1);

        // multiply by main gate selector(z)
        res.point_mul_assign(proof.gate_selector_values_at_z[0]); // these is only one explicitly opened selector

        uint256 current_alpha = 1;

        // calculate scalar contribution from the range check gate
        (current_alpha, tmp_fr) = add_contribution_from_range_constraint_gates(state, proof, current_alpha);

        tmp_g1 = vk.gate_selector_commitments[1].point_mul(tmp_fr); // selector commitment for range constraint gate * scalar
        res.point_add_assign(tmp_g1);

        // proceed as normal to copy permutation
        current_alpha = current_alpha.mulm(state.alpha); // alpha^5

        uint256 alpha_for_grand_product = current_alpha;

        // z * non_res * beta + gamma + a
        uint256 grand_product_part_at_z = state.z;
        grand_product_part_at_z = grand_product_part_at_z.mulm(state.beta);
        grand_product_part_at_z = grand_product_part_at_z.addm(proof.wire_values_at_z[0]);
        grand_product_part_at_z = grand_product_part_at_z.addm(state.gamma);
        for (uint256 i = 0; i < vk.copy_permutation_non_residues.length; i++) {
            tmp_fr = state.z;
            tmp_fr = tmp_fr.mulm(vk.copy_permutation_non_residues[i]);
            tmp_fr = tmp_fr.mulm(state.beta);
            tmp_fr = tmp_fr.addm(state.gamma);
            tmp_fr = tmp_fr.addm(proof.wire_values_at_z[i+1]);

            grand_product_part_at_z = grand_product_part_at_z.mulm(tmp_fr);
        }

        grand_product_part_at_z = grand_product_part_at_z.mulm(alpha_for_grand_product);

        // alpha^n & L_{0}(z), and we bump current_alpha
        current_alpha = current_alpha.mulm(state.alpha);

        tmp_fr = state.cached_lagrange_evals[0];
        tmp_fr = tmp_fr.mulm(current_alpha);

        grand_product_part_at_z = grand_product_part_at_z.addm(tmp_fr);

        // prefactor for grand_product(x) is complete

        // add to the linearization a part from the term
        // - (a(z) + beta*perm_a + gamma)*()*()*z(z*omega) * beta * perm_d(X)
        uint256 last_permutation_part_at_z = 1;
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            tmp_fr = state.beta;
            tmp_fr = tmp_fr.mulm(proof.permutation_polynomials_at_z[i]);
            tmp_fr = tmp_fr.addm(state.gamma);
            tmp_fr = tmp_fr.addm(proof.wire_values_at_z[i]);

            last_permutation_part_at_z = last_permutation_part_at_z.mulm(tmp_fr);
        }

        last_permutation_part_at_z = last_permutation_part_at_z.mulm(state.beta);
        last_permutation_part_at_z = last_permutation_part_at_z.mulm(proof.copy_grand_product_at_z_omega);
        last_permutation_part_at_z = last_permutation_part_at_z.mulm(alpha_for_grand_product); // we multiply by the power of alpha from the argument

        // actually multiply prefactors by z(x) and perm_d(x) and combine them
        tmp_g1 = proof.copy_permutation_grand_product_commitment.point_mul(grand_product_part_at_z);
        tmp_g1.point_sub_assign(vk.copy_permutation_commitments[STATE_WIDTH - 1].point_mul(last_permutation_part_at_z));

        res.point_add_assign(tmp_g1);
        // multiply them by v immedately as linearization has a factor of v^1
        res.point_mul_assign(state.v);
        // res now contains contribution from the gates linearization and
        // copy permutation part

        // now we need to add a part that is the rest
        // for z(x*omega):
        // - (a(z) + beta*perm_a + gamma)*()*()*(d(z) + gamma) * z(x*omega)
    }

    function aggregate_commitments(
        PartialVerifierState memory state,
        Proof memory proof,
        VerificationKey memory vk
    ) internal view returns (PairingsBn254.G1Point[2] memory res) {
        PairingsBn254.G1Point memory d = reconstruct_linearization_commitment(state, proof, vk);

        uint256 z_in_domain_size = state.z.pow(vk.domain_size);

        PairingsBn254.G1Point memory tmp_g1 = PairingsBn254.P1();

        uint256 aggregation_challenge = 1;

        PairingsBn254.G1Point memory commitment_aggregation = PairingsBn254.copy_g1(proof.quotient_poly_commitments[0]);
        uint256 tmp_fr = 1;
        for (uint i = 1; i < proof.quotient_poly_commitments.length; i++) {
            tmp_fr = tmp_fr.mulm(z_in_domain_size);
            tmp_g1 = proof.quotient_poly_commitments[i].point_mul(tmp_fr);
            commitment_aggregation.point_add_assign(tmp_g1);
        }

        aggregation_challenge = aggregation_challenge.mulm(state.v);
        commitment_aggregation.point_add_assign(d);

        for (uint i = 0; i < proof.wire_commitments.length; i++) {
            aggregation_challenge = aggregation_challenge.mulm(state.v);
            tmp_g1 = proof.wire_commitments[i].point_mul(aggregation_challenge);
            commitment_aggregation.point_add_assign(tmp_g1);
        }

        for (uint i = 0; i < NUM_GATE_SELECTORS_OPENED_EXPLICITLY; i++) {
            aggregation_challenge = aggregation_challenge.mulm(state.v);
            tmp_g1 = vk.gate_selector_commitments[0].point_mul(aggregation_challenge);
            commitment_aggregation.point_add_assign(tmp_g1);
        }

        for (uint i = 0; i < vk.copy_permutation_commitments.length - 1; i++) {
            aggregation_challenge = aggregation_challenge.mulm(state.v);
            tmp_g1 = vk.copy_permutation_commitments[i].point_mul(aggregation_challenge);
            commitment_aggregation.point_add_assign(tmp_g1);
        }

        aggregation_challenge = aggregation_challenge.mulm(state.v);
       // now do prefactor for grand_product(x*omega)
        tmp_fr = aggregation_challenge;
        tmp_fr = tmp_fr.mulm(state.u);
        commitment_aggregation.point_add_assign(proof.copy_permutation_grand_product_commitment.point_mul(tmp_fr));

        aggregation_challenge = aggregation_challenge.mulm(state.v);

        tmp_fr = aggregation_challenge;
        tmp_fr = tmp_fr.mulm(state.u);
        tmp_g1 = proof.wire_commitments[STATE_WIDTH - 1].point_mul(tmp_fr);
        commitment_aggregation.point_add_assign(tmp_g1);

        // collect opening values
        aggregation_challenge = 1;

        uint256 aggregated_value = (proof.quotient_polynomial_at_z);

        aggregation_challenge = aggregation_challenge.mulm(state.v);

        tmp_fr = proof.linearization_polynomial_at_z;
        tmp_fr = tmp_fr.mulm(aggregation_challenge);
        aggregated_value = aggregated_value.addm(tmp_fr);

        for (uint i = 0; i < proof.wire_values_at_z.length; i++) {
            aggregation_challenge = aggregation_challenge.mulm(state.v);

            tmp_fr = proof.wire_values_at_z[i];
            tmp_fr = tmp_fr.mulm(aggregation_challenge);
            aggregated_value = aggregated_value.addm(tmp_fr);
        }

        for (uint i = 0; i < proof.gate_selector_values_at_z.length; i++) {
            aggregation_challenge = aggregation_challenge.mulm(state.v);
            tmp_fr = proof.gate_selector_values_at_z[i];
            tmp_fr = tmp_fr.mulm(aggregation_challenge);
            aggregated_value = aggregated_value.addm(tmp_fr);
        }

        for (uint i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            aggregation_challenge = aggregation_challenge.mulm(state.v);

            tmp_fr = proof.permutation_polynomials_at_z[i];
            tmp_fr = tmp_fr.mulm(aggregation_challenge);
            aggregated_value = aggregated_value.addm(tmp_fr);
        }

        aggregation_challenge = aggregation_challenge.mulm(state.v);

        tmp_fr = proof.copy_grand_product_at_z_omega;
        tmp_fr = tmp_fr.mulm(aggregation_challenge);
        tmp_fr = tmp_fr.mulm(state.u);
        aggregated_value = aggregated_value.addm(tmp_fr);

        aggregation_challenge = aggregation_challenge.mulm(state.v);

        tmp_fr = proof.wire_values_at_z_omega[0];
        tmp_fr = tmp_fr.mulm(aggregation_challenge);
        tmp_fr = tmp_fr.mulm(state.u);
        aggregated_value = aggregated_value.addm(tmp_fr);

        commitment_aggregation.point_sub_assign(PairingsBn254.P1().point_mul(aggregated_value));

        PairingsBn254.G1Point memory pair_with_generator = commitment_aggregation;
        pair_with_generator.point_add_assign(proof.opening_at_z_proof.point_mul(state.z));

        tmp_fr = state.z;
        tmp_fr = tmp_fr.mulm(vk.omega);
        tmp_fr = tmp_fr.mulm(state.u);
        pair_with_generator.point_add_assign(proof.opening_at_z_omega_proof.point_mul(tmp_fr));

        PairingsBn254.G1Point memory pair_with_x = proof.opening_at_z_omega_proof.point_mul(state.u);
        pair_with_x.point_add_assign(proof.opening_at_z_proof);
        pair_with_x.negate();

        res[0] = pair_with_generator;
        res[1] = pair_with_x;

        return res;
    }

    function verify_initial(
        PartialVerifierState memory state,
        Proof memory proof,
        VerificationKey memory vk
    ) internal view returns (bool) {
        require(proof.input_values.length == vk.num_inputs, "J"); /* !ERROR: Incorrect number of input */
        require(vk.num_inputs >= 1, "XI"); /* !ERROR: Number of inputs must be >= 1 */
        TranscriptLibrary.Transcript memory transcript = TranscriptLibrary.new_transcript();
        for (uint256 i = 0; i < vk.num_inputs; i++) {
            transcript.update_with_u256(proof.input_values[i]);
        }

        for (uint256 i = 0; i < proof.wire_commitments.length; i++) {
            transcript.update_with_g1(proof.wire_commitments[i]);
        }

        state.beta = transcript.get_challenge();
        state.gamma = transcript.get_challenge();

        transcript.update_with_g1(proof.copy_permutation_grand_product_commitment);
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

        if (!valid) {
            return false;
        }

        transcript.update_with_fr(proof.quotient_polynomial_at_z);

        for (uint256 i = 0; i < proof.wire_values_at_z.length; i++) {
            transcript.update_with_fr(proof.wire_values_at_z[i]);
        }

        for (uint256 i = 0; i < proof.wire_values_at_z_omega.length; i++) {
            transcript.update_with_fr(proof.wire_values_at_z_omega[i]);
        }

        transcript.update_with_fr(proof.gate_selector_values_at_z[0]);

        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            transcript.update_with_fr(proof.permutation_polynomials_at_z[i]);
        }

        transcript.update_with_fr(proof.copy_grand_product_at_z_omega);
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

    function aggregate_for_verification(Proof memory proof, VerificationKey memory vk) internal view returns (bool valid, PairingsBn254.G1Point[2] memory part) {
        PartialVerifierState memory state;

        valid = verify_initial(state, proof, vk);

        if (!valid) {
            return (valid, part);
        }

        part = aggregate_commitments(state, proof, vk);

        (valid, part);
    }

    function verify(Proof memory proof, VerificationKey memory vk) internal view returns (bool) {
        (bool valid, PairingsBn254.G1Point[2] memory recursive_proof_part) = aggregate_for_verification(proof, vk);
        if (!valid) {
            return false;
        }

        valid = PairingsBn254.pairingProd2(recursive_proof_part[0], PairingsBn254.P2(), recursive_proof_part[1], vk.g2_x);

        return valid;
    }

    function verify_recursive(
        Proof memory proof,
        VerificationKey memory vk,
        uint256 recursive_vks_root,
        uint8 max_valid_index,
        uint8[] memory recursive_vks_indexes,
        uint256[] memory individual_vks_inputs,
        uint256[16] memory subproofs_limbs
    ) internal view returns (bool) {
        (uint256 recursive_input, PairingsBn254.G1Point[2] memory aggregated_g1s) = reconstruct_recursive_public_input(
            recursive_vks_root, max_valid_index, recursive_vks_indexes,
            individual_vks_inputs, subproofs_limbs
        );

        require(recursive_input == proof.input_values[0], "K"); /* !ERROR: Inconsistent input/proof */

        (bool valid, PairingsBn254.G1Point[2] memory recursive_proof_part) = aggregate_for_verification(proof, vk);
        if (!valid) {
            return false;
        }

        recursive_proof_part[0].point_add_assign(aggregated_g1s[0]);
        recursive_proof_part[1].point_add_assign(aggregated_g1s[1]);

        valid = PairingsBn254.pairingProd2(recursive_proof_part[0], PairingsBn254.P2(), recursive_proof_part[1], vk.g2_x);

        return valid;
    }

    function reconstruct_recursive_public_input(
        uint256 recursive_vks_root,
        uint8 max_valid_index,
        uint8[] memory recursive_vks_indexes,
        uint256[] memory individual_vks_inputs,
        uint256[16] memory subproofs_aggregated
    ) internal pure returns (uint256 recursive_input, PairingsBn254.G1Point[2] memory reconstructed_g1s) {

        uint256 length = recursive_vks_indexes.length;
        require(length == individual_vks_inputs.length, "SQ"); /* !ERROR: Arguments' lengths don't match */

        bytes memory concatenated = new bytes(32 + length*33 + 32*16);
        uint256 base;
        uint8 index;
        uint256 input;

        assembly {
            base := add(0x20, concatenated)
            mstore(base, recursive_vks_root)
        }

        for (uint256 i = 0; i < length; i++) {
            index = recursive_vks_indexes[i];
            require(index <= max_valid_index, "V7"); /* !ERROR: Invalid index */
            concatenated[32+i] = bytes1(index);
        }

        base += 32 + length;

        for (uint256 i = 0; i < length; i++) {
            input = individual_vks_inputs[i];
            require(input < PairingsBn254.r_mod, "I"); /* !ERROR: Invalid input */
            assembly {
                mstore(base, input)
            }
            base += 32;
        }

        for (uint256 i = 0; i < 16; i++) {
            input = subproofs_aggregated[i];
            require(input < MAX_LIMB_VALUE, "Y3"); /* !ERROR: Invalid subproof */
            assembly {
                mstore(base, input)
            }
            base += 32;
        }

        bytes32 commitment = sha256(concatenated);
        recursive_input = uint256(commitment) & RECURSIVE_CIRCUIT_INPUT_COMMITMENT_MASK;

        reconstructed_g1s[0] = PairingsBn254.new_g1_checked(
            subproofs_aggregated[0] + (subproofs_aggregated[1] << LIMB_WIDTH) + (subproofs_aggregated[2] << 2*LIMB_WIDTH) + (subproofs_aggregated[3] << 3*LIMB_WIDTH),
            subproofs_aggregated[4] + (subproofs_aggregated[5] << LIMB_WIDTH) + (subproofs_aggregated[6] << 2*LIMB_WIDTH) + (subproofs_aggregated[7] << 3*LIMB_WIDTH)
        );

        reconstructed_g1s[1] = PairingsBn254.new_g1_checked(
            subproofs_aggregated[8] + (subproofs_aggregated[9] << LIMB_WIDTH) + (subproofs_aggregated[10] << 2*LIMB_WIDTH) + (subproofs_aggregated[11] << 3*LIMB_WIDTH),
            subproofs_aggregated[12] + (subproofs_aggregated[13] << LIMB_WIDTH) + (subproofs_aggregated[14] << 2*LIMB_WIDTH) + (subproofs_aggregated[15] << 3*LIMB_WIDTH)
        );

        return (recursive_input, reconstructed_g1s);
    }
}

contract KeyedVerifier is Plonk4VerifierWithAccessToDNext {
    uint256 constant SERIALIZED_PROOF_LENGTH = 34;

    function get_verification_key() internal pure returns(VerificationKey memory vk) {
        vk.domain_size = 4194304;
        vk.num_inputs = 1;
        vk.omega = 0x18c95f1ae6514e11a1b30fd7923947c5ffcec5347f16e91b4dd654168326bede;
        vk.gate_setup_commitments[0] = PairingsBn254.new_g1(
            0x191bf5eb647f7404974c2d544e233d8031c11c89293068a29518c4f3b5a82b87,
            0x0289da987409f93ba5c56239b117a4fd3a73bd16c0ba62df82fd0dbd07ba2b04
        );
        vk.gate_setup_commitments[1] = PairingsBn254.new_g1(
            0x03edf444cead2f57dbb73e9a0e2e98985899bb0793d40f34a76df3967a11815f,
            0x0a5a682b95673c2f0eb8fbc4de0b2f8bddd237922e2a15bbf8edc7f19bb84086
        );
        vk.gate_setup_commitments[2] = PairingsBn254.new_g1(
            0x1fda99eb18dc0eddcda13b81796b3c289b2b4f6f72ebf306334dc52d01206319,
            0x0b380b280f6e74e6e47381bef76d00f951d742fbbd9814f87fb8b062ce52b07f
        );
        vk.gate_setup_commitments[3] = PairingsBn254.new_g1(
            0x133299193c2b6934a5056a8e154525334072af2f9979219bff66d40ebffa5c2a,
            0x212c48e85f3f03383beca295d1980e8d4aa8ef3e32f9a6cdefa300feae4aafd9
        );
        vk.gate_setup_commitments[4] = PairingsBn254.new_g1(
            0x058441552b7a4ca0f6ee97be2f26ff6ca4a4ea7f4d90096dd6661f42d40a1d9a,
            0x13ab4586de2cc5d1078a56758f6f5487bc1106d9f967600933871d6d7a81be83
        );
        vk.gate_setup_commitments[5] = PairingsBn254.new_g1(
            0x23a87dffd37ecbd5578a600effce1d22c79c9e869c3804420911be175c6fe665,
            0x0e8f6b12c32d24dcb05721a416a6ea229b6b22d2b2ffaaab308ebf873b2a86a9
        );
        vk.gate_setup_commitments[6] = PairingsBn254.new_g1(
            0x0e4d82b94eb786f366e56eea23570a83fbfde56fb4b5328dc65e0ba416a9d66e,
            0x23b6a34454e1ef57543f506e032fc3666583245fa8d0e1972d1ffc2a03b39bba
        );

        vk.gate_selector_commitments[0] = PairingsBn254.new_g1(
            0x033ced9a86f0149c5e11ccb5f83e390b730ce1d03352ea111d2ec19586c05a51,
            0x1e5c131310795ba1e6bbaf8c0700eb9a03929fcc39491757d596780de8147c08
        );
        vk.gate_selector_commitments[1] = PairingsBn254.new_g1(
            0x070cf4fdda67b7e084b2c96ac73fb7ff9119fd447c5ab7e5daabfa5dd5f3395e,
            0x06047d1a934d58f61f345b91e9aebab1092102ebcab0ba9196119784a881767d
        );

        vk.copy_permutation_commitments[0] = PairingsBn254.new_g1(
            0x0b1b50cee2777e03088ffaffcae0d77136274d721eb7fe5ef97ca72e6ec235f5,
            0x145cb123ae891ae63f525d447d52f90206824d3616ccc85ebb5eff8f2ff109b2
        );
        vk.copy_permutation_commitments[1] = PairingsBn254.new_g1(
            0x0374b6a602699e56ac5ebf8f7b1ee8f432a54d8f133ff0d0b4c7a68efc9e9513,
            0x24f22112315c574096f642b60b372494a661999db5b1b57c2a0af1f93329f548
        );
        vk.copy_permutation_commitments[2] = PairingsBn254.new_g1(
            0x2663bb35e725e7e6b473b753032ae175193674388311e82a2bd1c2a008f0dd40,
            0x02965c14a41c8846d4e0abe0de3786209aaa4db2689a9401fb0f1a562afca435
        );
        vk.copy_permutation_commitments[3] = PairingsBn254.new_g1(
            0x2dd86ba3b5fceea533ae2feba3f01e97261faf9e458501eafa5dd608a50049da,
            0x26e06bb657495b8d898ca7065ebc2f47f9f53724ff684f0b9d633e0c69c523ed
        );

        vk.copy_permutation_non_residues[0] =
            0x0000000000000000000000000000000000000000000000000000000000000005;
        vk.copy_permutation_non_residues[1] =
            0x0000000000000000000000000000000000000000000000000000000000000007;
        vk.copy_permutation_non_residues[2] =
            0x000000000000000000000000000000000000000000000000000000000000000a;

        vk.g2_x = PairingsBn254.new_g2(
            [0x260e01b251f6f1c7e7ff4e580791dee8ea51d87a358e038b4efe30fac09383c1,
             0x0118c4d5b837bcc2bc89b5b398b5974e9f5944073b32078b7e231fec938883b0],
            [0x04fc6369f7110fe3d25156c1bb9a72859cf2a04641f99ba4ee413c80da6a5fe4,
             0x22febda3c0c0632a56475b4214e5615e11e6dd3f96e6cea2854a87d4dacc5e55]
        );
    }

    function deserialize_proof(
        uint256[] memory public_inputs,
        uint256[] memory serialized_proof
    ) internal pure returns(Proof memory proof) {
        require(serialized_proof.length == SERIALIZED_PROOF_LENGTH, "2J"); /* !ERROR: Invalid proof length */
        proof.input_values = new uint256[](public_inputs.length);
        for (uint256 i = 0; i < public_inputs.length; i++) {
            proof.input_values[i] = public_inputs[i];
        }

        uint256 j = 0;
        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.wire_commitments[i] = PairingsBn254.new_g1_checked(
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

        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.quotient_poly_commitments[i] = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
            );

            j += 2;
        }

        for (uint256 i = 0; i < STATE_WIDTH; i++) {
            proof.wire_values_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }

        for (uint256 i = 0; i < proof.wire_values_at_z_omega.length; i++) {
            proof.wire_values_at_z_omega[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }

        for (uint256 i = 0; i < proof.gate_selector_values_at_z.length; i++) {
            proof.gate_selector_values_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }

        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            proof.permutation_polynomials_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j += 1;
        }


        proof.copy_grand_product_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j += 1;

        proof.quotient_polynomial_at_z = PairingsBn254.new_fr(
            serialized_proof[j]
        );

        j += 1;

        proof.linearization_polynomial_at_z = PairingsBn254.new_fr(
            serialized_proof[j]
        );

        j += 1;

        proof.opening_at_z_proof = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;

        proof.opening_at_z_omega_proof = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j+1]
        );
    }

    function verify_serialized_proof(
        uint256[] memory public_inputs,
        uint256[] memory serialized_proof
    ) public view returns (bool) {
        VerificationKey memory vk = get_verification_key();
        require(vk.num_inputs == public_inputs.length, "G"); /* !ERROR: Invalid number of inputs */

        Proof memory proof = deserialize_proof(public_inputs, serialized_proof);

        bool valid = verify(proof, vk);

        return valid;
    }

    function verify_serialized_proof_with_recursion(
        uint256[] memory public_inputs,
        uint256[] memory serialized_proof,
        uint256 recursive_vks_root,
        uint8 max_valid_index,
        uint8[] memory recursive_vks_indexes,
        uint256[] memory individual_vks_inputs,
        uint256[16] memory subproofs_limbs
    ) public view returns (bool) {
        VerificationKey memory vk = get_verification_key();
        require(vk.num_inputs == public_inputs.length, "R"); /* !ERROR: Invalid length of input */

        Proof memory proof = deserialize_proof(public_inputs, serialized_proof);

        bool valid = verify_recursive(proof, vk, recursive_vks_root, max_valid_index, recursive_vks_indexes, individual_vks_inputs, subproofs_limbs);

        return valid;
    }
}
