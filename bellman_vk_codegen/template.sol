pragma solidity >=0.6.0 <0.7.0;

library Pairings {
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
        uint256 r = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
        require(fr < r);
        return Fr({value: fr});
    }
    
    function copy(Fr memory self) internal pure returns (Fr memory n) {
        n.value = self.value;
    }
    
    function assign(Fr memory self, Fr memory other) internal pure {
        self.value = other.value;
    }
    
    function inverse(Fr memory fr) internal view returns (Fr memory) {
        return pow(fr, r_mod-2);
    }
    
    function add_assign(Fr memory self, Fr memory other) internal pure {
        uint256 r = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
        self.value = addmod(self.value, other.value, r);
    }
    
    function sub_assign(Fr memory self, Fr memory other) internal pure {
        uint256 r = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
        self.value = addmod(self.value, r - other.value, r);
    }
    
    function mul_assign(Fr memory self, Fr memory other) internal pure {
        uint256 r = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
        self.value = mulmod(self.value, other.value, r);
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
        uint256 q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
        if (self.X == 0 && self.Y == 0)
            return;
        self.Y = q - self.Y;
    }
    
    function is_infinity(G1Point memory p) internal pure returns (bool) {
        if (p.X == 0) {
            if (p.Y == 0 || p.Y == 1) {
                return true;
            }
        }
        return false;
    }
    
    function normalize(G1Point memory self) internal pure {
        if (self.X == 0) {
            if (self.Y == 1) {
                self.Y = 0;
            }
        }
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

        uint256 q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
        if (p1.X == 0 && p1.Y == 0) {
            p1.X = p2.X;
            p1.Y = p2.Y;
            return;
        } else if (p2.X == 0 && p2.Y == 0) {
            return;
        } else {
            input[0] = p1.X;
            input[1] = p1.Y;
            input[2] = p2.X;
            input[3] = q - p2.Y;
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
    
    function multiexp(G1Point[] memory bases, Fr[] memory values)
        internal view returns (G1Point memory r)
    {
        require(bases.length == values.length);
        require(bases.length >= 2);
        
        r = pointMul(bases[0], values[0]);
        G1Point memory tmp = pointMul(bases[1], values[1]);
        r = pointAdd(r, tmp);
        
        for (uint i = 2; i < bases.length; i++) {
            tmp = pointMul(bases[i], values[i]);
            r = pointAdd(r, tmp);
        }
        
        return r;
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
    uint256 constant mask = 0xe000000000000000000000000000000000000000000000000000000000000000;
    
    struct Transcript {
        bytes32 state;
    }
    
    function update_with_fr(Transcript memory self, Pairings.Fr memory value) internal pure {
        bytes32 new_state = keccak256(abi.encodePacked(self.state, value.value));
        self.state = new_state;
    }
    
    function update_with_g1(Transcript memory self, Pairings.G1Point memory p) internal pure {
        bytes32 new_state = keccak256(abi.encodePacked(self.state, p.X));
        new_state = keccak256(abi.encodePacked(new_state, p.Y));
        self.state = new_state;
    }
    
    function update_with_u256(Transcript memory self, uint256 value) internal pure {
        bytes32 new_state = keccak256(abi.encodePacked(self.state, value));
        self.state = new_state;
    }
    
    function get_challenge(Transcript memory self) internal pure returns(Pairings.Fr memory challenge) {
        challenge = Pairings.Fr({value: uint256(self.state) & (~mask)});
        bytes32 new_state = keccak256(abi.encodePacked(self.state));
        self.state = new_state;
    }
}

contract Plonk4Verifier {
    using Pairings for Pairings.G1Point;
    using Pairings for Pairings.G2Point;
    using Pairings for Pairings.Fr;
    
    using TranscriptLibrary for TranscriptLibrary.Transcript;
    
    event DebugU256(uint256);
    event DebugBytes32(bytes32);
    
    struct VerificationKey {
        uint256 domain_size;
        uint256 num_inputs;
        Pairings.Fr omega;
        Pairings.G1Point[] selector_commitments;
        Pairings.G1Point[] next_step_selector_commitments;
        Pairings.G1Point[] permutation_commitments;
        Pairings.Fr[] permutation_non_residues;
        Pairings.G2Point g2_x;
    }
    
    struct Proof {
        uint256[] input_values;
        Pairings.G1Point[] wire_commitments;
        Pairings.G1Point grand_product_commitment;
        Pairings.G1Point[] quotient_poly_commitments;
        Pairings.Fr[] wire_values_at_z;
        Pairings.Fr[] wire_values_at_z_omega;
        Pairings.Fr grand_product_at_z_omega;
        Pairings.Fr quotient_polynomial_at_z;
        Pairings.Fr linearization_polynomial_at_z;
        Pairings.Fr[] permutation_polynomials_at_z;
    
        Pairings.G1Point opening_at_z_proof;
        Pairings.G1Point opening_at_z_omega_proof;
    }
    
    struct PartialVerifierState {
        Pairings.Fr alpha;
        Pairings.Fr beta;
        Pairings.Fr gamma;
        Pairings.Fr v;
        Pairings.Fr u;
        Pairings.Fr z;
        Pairings.Fr[] cached_lagrange_evals;
    }
    
    uint256 constant state_width = 4;
    
    function get_verification_key() internal pure returns(VerificationKey memory vk) {
        vk.domain_size = {{domain_size}};
        vk.num_inputs = {{num_inputs}};
        vk.omega = Pairings.newFr({{omega}});
        vk.selector_commitments = new Pairings.G1Point[](state_width+2);
        vk.selector_commitments[0] = Pairings.newG1(
            {{selector_commitment_0_0}},
            {{selector_commitment_0_1}}
        );
        vk.selector_commitments[1] = Pairings.newG1(
            {{selector_commitment_1_0}},
            {{selector_commitment_1_1}}
        );
        vk.selector_commitments[2] = Pairings.newG1(
            {{selector_commitment_2_0}},
            {{selector_commitment_2_1}}
        );
        vk.selector_commitments[3] = Pairings.newG1(
            {{selector_commitment_3_0}},
            {{selector_commitment_3_1}}
        );
        vk.selector_commitments[4] = Pairings.newG1(
            {{selector_commitment_4_0}},
            {{selector_commitment_4_1}}
        );
        vk.selector_commitments[5] = Pairings.newG1(
            {{selector_commitment_5_0}},
            {{selector_commitment_5_1}}
        );
        
        vk.next_step_selector_commitments = new Pairings.G1Point[](1);
        vk.next_step_selector_commitments[0] = Pairings.newG1(
            {{next_step_selector_commitment_0_0}},
            {{next_step_selector_commitment_0_1}}
        );
        
        vk.permutation_commitments = new Pairings.G1Point[](state_width);
         vk.permutation_commitments[0] = Pairings.newG1(
            {{permutation_commitment_0_0}},
            {{permutation_commitment_0_1}}
        );
        vk.permutation_commitments[1] = Pairings.newG1(
            {{permutation_commitment_1_0}},
            {{permutation_commitment_1_1}}
        );
        vk.permutation_commitments[2] = Pairings.newG1(
            {{permutation_commitment_2_0}},
            {{permutation_commitment_2_1}}
        );
        vk.permutation_commitments[3] = Pairings.newG1(
            {{permutation_commitment_3_0}},
            {{permutation_commitment_3_1}}
        );
        
        vk.permutation_non_residues = new Pairings.Fr[](state_width-1);
        vk.permutation_non_residues[0] = Pairings.newFr(
            {{permutation_non_residue_0}}
        );
        vk.permutation_non_residues[1] = Pairings.newFr(
            {{permutation_non_residue_1}}
        );
        vk.permutation_non_residues[2] = Pairings.newFr(
            {{permutation_non_residue_2}}
        );
        
        vk.g2_x = Pairings.newG2(
            [{{g2_x_x_c1}},
             {{g2_x_x_c0}}],
            [{{g2_x_y_c1}},
             {{g2_x_y_c0}}]
        );
    }
    
    function get_dummy_proof() internal pure returns(Proof memory proof) {
        proof.input_values = new uint256[](1);
        proof.input_values[0] = 1;
        proof.wire_commitments = new Pairings.G1Point[](state_width);
        proof.wire_commitments[0] = Pairings.newG1(
            0x1b5e58c31beca50e963688ba8d6be08b1c90e964135450d882876c8a5e8bc613,
            0x2d2321db916c233e56596bba357913cdad1e962e9264f9c74773af210a5ecd29
        );
        proof.wire_commitments[1] = Pairings.newG1(
            0x2c4951b832621f1b7315181f70905ba8cd4865753d98da5c28a474092edadbda,
            0x27b14de154fd0af1394e5ca8f8c64032ee91cd9e046cbd7a24aeed18db9c3334
        );
        proof.wire_commitments[2] = Pairings.newG1(
            0x1a491b60cdae811dc03c00655a64073cf536115324b340895c6db7e387eea681,
            0x06e63b071abae5ad0dc281f103e234f34535d795699dcc6a895d5964422a5952
        );
        proof.wire_commitments[3] = Pairings.newG1(
            0x0000000000000000000000000000000000000000000000000000000000000000,
            0x0000000000000000000000000000000000000000000000000000000000000000
        );
        
        proof.grand_product_commitment = Pairings.newG1(
            0x161b6ef60ada7b486c9bfbc23350365422fc3eff69aa0716c1319b696f41dc0b,
            0x1cbc0c1bd268e96ac40a377067ee2ca383983870ce555db61ec2ab25fa3996c4
        );
        
        proof.quotient_poly_commitments = new Pairings.G1Point[](state_width);
        proof.quotient_poly_commitments[0] = Pairings.newG1(
            0x1d8f5d10f66ea913f4abd01d178e6c0367f6ed7e6ad86f0637eec170043b794c,
            0x16339d848618af3adfc06827fe67a3b7c65a989b442d6e4266056005052bb72c
        );
        proof.quotient_poly_commitments[1] = Pairings.newG1(
            0x2012ea21d65fe7bd0be88491eba50e779becced0c72030babce4ea652fd5e92a,
            0x28c63d4ef6bbceb1602d5e4df8c7fe8e33a2360241c81fc320f16910bb05ed09
        );
        proof.quotient_poly_commitments[2] = Pairings.newG1(
            0x0041c38c5baae54d3e9cdd56a85e88ed976f3b082bcfef2067658886372067ba,
            0x062c2b76fa36f2bf0bcdaa17a27e8c10053bec7870c599086a6338d01dbd3262
        );
        proof.quotient_poly_commitments[3] = Pairings.newG1(
            0x0000000000000000000000000000000000000000000000000000000000000000,
            0x0000000000000000000000000000000000000000000000000000000000000000
        );
        
        proof.wire_values_at_z = new Pairings.Fr[](state_width);
        proof.wire_values_at_z[0] = Pairings.newFr(
            0x1b6902ca32bc26a34f6b846eeb93ac69fa4581617da3b318460eca73dd2d18e0
        );
        proof.wire_values_at_z[1] = Pairings.newFr(
            0x2c9adab6d77eb011f3e3e511f2f0df1e8fc695601c7775f8dfa1d881b244dc4c
        );
        proof.wire_values_at_z[2] = Pairings.newFr(
            0x2a139912dcbb852a87c257a100793c63da16ca0fd146177e4f47fbd066706f2d
        );
        proof.wire_values_at_z[3] = Pairings.newFr(
            0x0000000000000000000000000000000000000000000000000000000000000000
        );
        
        proof.wire_values_at_z_omega = new Pairings.Fr[](1);
        proof.wire_values_at_z_omega[0] = Pairings.newFr(
            0x0000000000000000000000000000000000000000000000000000000000000000
        );
        
        proof.grand_product_at_z_omega = Pairings.newFr(
            0x09334ae6db6affb8c1efb113552128645b728f9b64a4b28466243a8fdf31dcc3
        );
        
        proof.quotient_polynomial_at_z = Pairings.newFr(
            0x06dcd12571bf33be988cfb99610eb0c3e6fd135cee07799afb177a9a710ac7b3
        );
        
        proof.linearization_polynomial_at_z = Pairings.newFr(
            0x013c7a8c3d10661074f62207bacc0d98105182b29a9b6f077c7561cf4aea90f7
        );
        
        proof.permutation_polynomials_at_z = new Pairings.Fr[](state_width - 1);
        proof.permutation_polynomials_at_z[0] = Pairings.newFr(
            0x1d8cdc7c9e8e4c7b07d69f62746da251c5de4f0e677f4b157621e47a097c0f97
        );
        proof.permutation_polynomials_at_z[1] = Pairings.newFr(
            0x284d14ca9b828c96398366151333fac8ca4af048e7a5bf424d5c77a2731a8132
        );
        proof.permutation_polynomials_at_z[2] = Pairings.newFr(
            0x00f766126ed12460d91886cccd2e5f467f1856aec6a8f39c3b732dd4883548ba
        );
        
        proof.opening_at_z_proof = Pairings.newG1(
            0x0f0ec46b791690bfcc86e7e98c65b431e9ff623bed4ba78ef5ba1ff9c9e24468,
            0x23b2037afd719f15f230aabc2631d8dca79170d013cd0182608e36fb4f18da0f
        );
        proof.opening_at_z_omega_proof = Pairings.newG1(
            0x14cdb2e64f4deef40be7a70381246a752698e8b798d3f00f4dee0e2e216acc1d,
            0x199a18224eb1774d3646008365ec59c333b385262d5a32e9131cdd1493ddeb7f
        );
    }
    
    function evaluate_lagrange_poly(
        uint256 poly_num, 
        uint256 domain_size, 
        Pairings.Fr memory omega, 
        Pairings.Fr memory at
    ) internal view returns (Pairings.Fr memory res) {
        require(poly_num < domain_size);
        Pairings.Fr memory one = Pairings.newFr(1);
        Pairings.Fr memory omega_power = omega.pow(poly_num);
        res = at.pow(domain_size);
        res.sub_assign(one);
        res.mul_assign(omega_power);
        
        Pairings.Fr memory den = Pairings.copy(at);
        den.sub_assign(omega_power);
        den.mul_assign(Pairings.newFr(domain_size));
        
        den = den.inverse();
        
        res.mul_assign(den);
    }
    
    function batch_evaluate_lagrange_poly(
        uint256[] memory poly_nums, 
        uint256 domain_size, 
        Pairings.Fr memory omega, 
        Pairings.Fr memory at
    ) internal view returns (Pairings.Fr[] memory res) {
        Pairings.Fr memory one = Pairings.newFr(1);
        Pairings.Fr memory tmp_1 = Pairings.newFr(0);
        Pairings.Fr memory tmp_2 = Pairings.newFr(domain_size);
        Pairings.Fr memory power_of_z = at.pow(domain_size);
        Pairings.Fr[] memory nums = new Pairings.Fr[](poly_nums.length);
        Pairings.Fr[] memory dens = new Pairings.Fr[](poly_nums.length);
        // numerators in a form omega^i * (z^n - 1)
        // denoms in a form (z - omega^i) * N
        for (uint i = 0; i < poly_nums.length; i++) {
            tmp_1 = omega.pow(poly_nums[i]); // power of omega
            nums[i].assign(power_of_z);
            nums[i].sub_assign(one);
            nums[i].mul_assign(tmp_1);
            
            dens[i].assign(at); // (X - omega^i) * N
            dens[i].sub_assign(tmp_1); 
            dens[i].mul_assign(tmp_2); // mul by domain size
        }
        
        Pairings.Fr[] memory partial_products = new Pairings.Fr[](poly_nums.length);
        partial_products[0].assign(Pairings.newFr(1));
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
        Pairings.Fr memory at
    ) internal view returns (Pairings.Fr memory res) {
        res = at.pow(domain_size);
        res.sub_assign(Pairings.newFr(1));
    }
    
    function verify_at_z(
        PartialVerifierState memory state,
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (bool) {
        Pairings.Fr memory lhs = evaluate_vanishing(vk.domain_size, state.z);
        lhs.mul_assign(proof.quotient_polynomial_at_z);
    
        Pairings.Fr memory quotient_challenge = Pairings.newFr(1);
        Pairings.Fr memory rhs = Pairings.copy(proof.linearization_polynomial_at_z);
        
        // public inputs
        Pairings.Fr memory tmp = Pairings.newFr(0);
        for (uint256 i = 0; i < proof.input_values.length; i++) {
            tmp.assign(state.cached_lagrange_evals[i]);
            tmp.mul_assign(Pairings.newFr(proof.input_values[i]));
            rhs.add_assign(tmp);
        }
        
        quotient_challenge.mul_assign(state.alpha);
        
        Pairings.Fr memory z_part = Pairings.copy(proof.grand_product_at_z_omega);
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            tmp.assign(proof.permutation_polynomials_at_z[i]);
            tmp.mul_assign(state.beta);
            tmp.add_assign(state.gamma);
            tmp.add_assign(proof.wire_values_at_z[i]);
            
            z_part.mul_assign(tmp);
        }
        
        tmp.assign(state.gamma);
        tmp.add_assign(proof.wire_values_at_z[proof.permutation_polynomials_at_z.length]);
        
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
    ) internal view returns (Pairings.G1Point memory res) {
        uint256 state_w = state_width;
        uint256 power_for_z_omega_opening = 1 + 1 + state_w + state_w - 1;
        res = Pairings.copyG1(vk.selector_commitments[state_w + 1]);
                
        Pairings.G1Point memory tmp_g1 = Pairings.P1();
        Pairings.Fr memory tmp_fr = Pairings.newFr(0);
        
        // addition gates
        for (uint256 i = 0; i < state_w; i++) {
            tmp_g1 = vk.selector_commitments[i].pointMul(proof.wire_values_at_z[i]);
            res.point_add_assign(tmp_g1);
        }
        
        // multiplication gate
        tmp_fr.assign(proof.wire_values_at_z[0]);
        tmp_fr.mul_assign(proof.wire_values_at_z[1]);
        tmp_g1 = vk.selector_commitments[state_w].pointMul(tmp_fr);
        res.point_add_assign(tmp_g1);
        
        // d_next
        tmp_g1 = vk.next_step_selector_commitments[0].pointMul(proof.wire_values_at_z_omega[0]);
        res.point_add_assign(tmp_g1);
        
        // z * non_res * beta + gamma + a
        Pairings.Fr memory grand_product_part_at_z = Pairings.copy(state.z);
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
        
        Pairings.Fr memory grand_product_part_at_z_omega = state.v.pow(power_for_z_omega_opening);
        grand_product_part_at_z_omega.mul_assign(state.u);
        
        Pairings.Fr memory last_permutation_part_at_z = Pairings.newFr(1);
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
        tmp_g1.point_sub_assign(vk.permutation_commitments[state_width - 1].pointMul(last_permutation_part_at_z));

        res.point_add_assign(tmp_g1);
        res.point_mul_assign(state.v);
        
        res.point_add_assign(proof.grand_product_commitment.pointMul(grand_product_part_at_z_omega));
    }
    
    function verify_commitments(
        PartialVerifierState memory state,
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (bool) {
        Pairings.G1Point memory d = reconstruct_d(state, proof, vk);
        
        Pairings.Fr memory z_in_domain_size = state.z.pow(vk.domain_size);
        
        Pairings.G1Point memory tmp_g1 = Pairings.P1();

        Pairings.Fr memory aggregation_challenge = Pairings.newFr(1);
        
        Pairings.G1Point memory commitment_aggregation = Pairings.copyG1(proof.quotient_poly_commitments[0]);
        Pairings.Fr memory tmp_fr = Pairings.newFr(1);
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
        tmp_g1 = proof.wire_commitments[state_width - 1].pointMul(tmp_fr);
        commitment_aggregation.point_add_assign(tmp_g1);
        
        // collect opening values
        aggregation_challenge = Pairings.newFr(1);
        
        Pairings.Fr memory aggregated_value = Pairings.copy(proof.quotient_polynomial_at_z);
        
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
        
        commitment_aggregation.point_sub_assign(Pairings.P1().pointMul(aggregated_value));
        
        Pairings.G1Point memory pair_with_generator = commitment_aggregation;
        pair_with_generator.point_add_assign(proof.opening_at_z_proof.pointMul(state.z));
        
        tmp_fr.assign(state.z);
        tmp_fr.mul_assign(vk.omega);
        tmp_fr.mul_assign(state.u);
        pair_with_generator.point_add_assign(proof.opening_at_z_omega_proof.pointMul(tmp_fr));
        
        Pairings.G1Point memory pair_with_x = proof.opening_at_z_omega_proof.pointMul(state.u);
        pair_with_x.point_add_assign(proof.opening_at_z_proof);
        pair_with_x.negate();
        
        return Pairings.pairingProd2(pair_with_generator, Pairings.P2(), pair_with_x, vk.g2_x);
    }
    
    function verify_initial(
        PartialVerifierState memory state, 
        Proof memory proof, 
        VerificationKey memory vk
    ) internal view returns (bool) {
        require(proof.input_values.length == vk.num_inputs);
        require(vk.num_inputs >= 1);
        TranscriptLibrary.Transcript memory transcript = TranscriptLibrary.Transcript({state: keccak256(abi.encodePacked(proof.input_values[0]))});
        for (uint256 i = 1; i < vk.num_inputs; i++) {
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
        
        state.cached_lagrange_evals = batch_evaluate_lagrange_poly(
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
    
    function verify(Proof memory proof, VerificationKey memory vk) internal view returns (bool) {
        PartialVerifierState memory state;
        
        bool valid = verify_initial(state, proof, vk);
        
        if (valid == false) {
            return false;
        }
        
        valid = verify_commitments(state, proof, vk);
        
        return valid;
    }

    function deserialize_proof(uint256 expected_inputs, uint256[] memory public_inputs, uint256[] memory serialized_proof) internal pure returns(Proof memory proof) {
        assert(expected_inputs == public_inputs.length);
        proof.input_values = new uint256[](expected_inputs);
        for (uint256 i = 0; i < expected_inputs; i++) {
            proof.input_values[i] = public_inputs[i];
        }
 
        uint256 j = 0;
        proof.wire_commitments = new Pairings.G1Point[](state_width);
        for (uint256 i = 0; i < state_width; i++) {
            proof.wire_commitments[i] = Pairings.newG1(
                serialized_proof[j],
                serialized_proof[j+1]
            );

            j += 2;
        }
        
        proof.grand_product_commitment = Pairings.newG1(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;
        
        proof.quotient_poly_commitments = new Pairings.G1Point[](state_width);
        for (uint256 i = 0; i < state_width; i++) {
            proof.quotient_poly_commitments[i] = Pairings.newG1(
                serialized_proof[j],
                serialized_proof[j+1]
            );

            j += 2;
        }
        
        proof.wire_values_at_z = new Pairings.Fr[](state_width);
        for (uint256 i = 0; i < state_width; i++) {
            proof.wire_values_at_z[i] = Pairings.newFr(
                serialized_proof[j]
            );

            j += 1;
        }
        
        proof.wire_values_at_z_omega = new Pairings.Fr[](1);
        for (uint256 i = 0; i < proof.wire_values_at_z_omega.length; i++) {
            proof.wire_values_at_z_omega[i] = Pairings.newFr(
                serialized_proof[j]
            );

            j += 1;
        }
        
        proof.grand_product_at_z_omega = Pairings.newFr(
                serialized_proof[j]
            );

        j += 1;

        proof.quotient_polynomial_at_z = Pairings.newFr(
            serialized_proof[j]
        );

        j += 1;

        proof.linearization_polynomial_at_z = Pairings.newFr(
            serialized_proof[j]
        );

        j += 1;
    
        proof.permutation_polynomials_at_z = new Pairings.Fr[](state_width - 1);
        for (uint256 i = 0; i < proof.permutation_polynomials_at_z.length; i++) {
            proof.permutation_polynomials_at_z[i] = Pairings.newFr(
                serialized_proof[j]
            );

            j += 1;
        }

        proof.opening_at_z_proof = Pairings.newG1(
                serialized_proof[j],
                serialized_proof[j+1]
        );
        j += 2;

        proof.opening_at_z_omega_proof = Pairings.newG1(
                serialized_proof[j],
                serialized_proof[j+1]
        );
    
        j += 2;
        assert(j == serialized_proof.length);
    }
    
    function test() public returns (bool) {
        bool valid = verify(get_dummy_proof(), get_verification_key());
        if (valid == false) {
            emit DebugBytes32(bytes32(uint256(0)));
        } else {
            emit DebugBytes32(bytes32(uint256(1)));
        }
    }

    function verity(
        uint256[] memory public_inputs, 
        uint256[] memory serialized_proof
    ) public view returns (bool) {
        VerificationKey memory vk = get_verification_key();
        uint256 expected_inputs = vk.num_inputs;

        Proof memory proof = deserialize_proof(expected_inputs, public_inputs, serialized_proof);

        bool valid = verify(proof, vk);

        return valid;
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
}
