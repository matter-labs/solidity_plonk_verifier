import { BigNumber } from "ethers";
import {readFileSync} from "fs";
import {utils} from "ethers"
const FR_SIZE: number = 32; // 

class Point{
    x: BigNumber;
    y: BigNumber;

    constructor(x: BigNumber, y: BigNumber){
        this.x = x;
        this.y = y;
    }

    static zero(): Point {
        return new Point(BigNumber.from(0), BigNumber.from(0));
    }

    encode(): BigNumber[]{
        return [this.x, this.y];
    }
}

export class Proof{
    n: BigNumber = BigNumber.from(0);
    inputs: BigNumber[] = [];
    state_polys_commitments: Point[] = [];
    witness_polys_commitments: Point[] = [];
    copy_permutation_grand_product_commitment: Point = Point.zero();
    
    lookup_s_poly_commitment?: Point;
    lookup_grand_product_commitment?: Point;
    
    quotient_poly_parts_commitments: Point[] = [];

    state_polys_openings_at_z: BigNumber[] = [];
    state_polys_openings_at_dilations: [BigNumber, BigNumber, BigNumber][] = [];
    witness_polys_openings_at_z: BigNumber[] = [];
    witness_polys_openings_at_dilations: [BigNumber, BigNumber, BigNumber][] = [];
    
    gate_setup_openings_at_z: [BigNumber, BigNumber, BigNumber][] = [];
    gate_selectors_openings_at_z: [BigNumber, BigNumber][] = [];
    
    copy_permutation_polys_openings_at_z: BigNumber[] = [];
    copy_permutation_grand_product_opening_at_z_omega: BigNumber =  BigNumber.from(0);

    lookup_s_poly_opening_at_z_omega?: BigNumber;
    lookup_grand_product_opening_at_z_omega?: BigNumber;
    lookup_t_poly_opening_at_z?: BigNumber ;
    lookup_t_poly_opening_at_z_omega?: BigNumber;

    lookup_selector_poly_opening_at_z?: BigNumber;
    lookup_table_type_poly_opening_at_z?: BigNumber;

    quotient_poly_opening_at_z: BigNumber = BigNumber.from(0);
    linearization_poly_opening_at_z: BigNumber = BigNumber.from(0);

    opening_proof_at_z: Point = Point.zero();
    opening_proof_at_z_omega: Point = Point.zero();

    static read_from_file(filePath: string): Proof{
        let buffer = readFileSync(filePath);
        let buf = new Buf(buffer); 
        
        return this.read(buf);
    }

    static read(buf: Buf): Proof{
        let proof = new Proof();
        proof.n = buf.read_u64();
        proof.inputs = buf.read_fr_vec(); 
        proof.state_polys_commitments = buf.read_curve_affine_vector(); 
        proof.witness_polys_commitments = buf.read_curve_affine_vector(); 
        proof.copy_permutation_grand_product_commitment = buf.read_curve_affine(); 
        
        proof.lookup_s_poly_commitment = buf.read_optional_curve_affine(); 
        proof.lookup_grand_product_commitment = buf.read_optional_curve_affine(); 
        proof.quotient_poly_parts_commitments = buf.read_curve_affine_vector(); 

        proof.state_polys_openings_at_z = buf.read_fr_vec();
        proof.state_polys_openings_at_dilations = buf.read_tuple_with_two_indexes_vec();
        proof.witness_polys_openings_at_z = buf.read_fr_vec();
        proof.witness_polys_openings_at_dilations = buf.read_tuple_with_two_indexes_vec();
        proof.gate_setup_openings_at_z = buf.read_tuple_with_two_indexes_vec();
        proof.gate_selectors_openings_at_z = buf.read_tuple_with_one_index_vec();
        proof.copy_permutation_polys_openings_at_z = buf.read_fr_vec();
        proof.copy_permutation_grand_product_opening_at_z_omega = buf.read_fr();
        
        proof.lookup_s_poly_opening_at_z_omega = buf.read_optional_fr();
        proof.lookup_grand_product_opening_at_z_omega = buf.read_optional_fr();
        proof.lookup_t_poly_opening_at_z = buf.read_optional_fr();
        proof.lookup_t_poly_opening_at_z_omega = buf.read_optional_fr();
        proof.lookup_selector_poly_opening_at_z = buf.read_optional_fr();
        proof.lookup_table_type_poly_opening_at_z = buf.read_optional_fr();
        
        proof.quotient_poly_opening_at_z = buf.read_fr();
        proof.linearization_poly_opening_at_z = buf.read_fr();

        proof.opening_proof_at_z = buf.read_curve_affine();
        proof.opening_proof_at_z_omega = buf.read_curve_affine();

        return proof;
    }

    serialize(): [BigNumber[], BigNumber[]]{
        let encoding : BigNumber[] = [];      
        for(let el of this.state_polys_commitments){
            encoding.push(...el.encode());
        }
        for(let el of this.witness_polys_commitments){
            encoding.push(...el.encode());
        }
        encoding.push(...this.copy_permutation_grand_product_commitment.encode())
        if(this.lookup_s_poly_commitment){
            encoding.push(...this.lookup_s_poly_commitment.encode())
        }
        if(this.lookup_grand_product_commitment){
            encoding.push(...this.lookup_grand_product_commitment.encode())
        }
        if(this.quotient_poly_parts_commitments){
            for(let el of this.quotient_poly_parts_commitments){
                encoding.push(...el.encode())
            }            
        }
        encoding.push(...this.state_polys_openings_at_z);
        for(let el of this.state_polys_openings_at_dilations){
            encoding.push(el[2])
        }
        encoding.push(...this.witness_polys_openings_at_z);
        for(let el of this.witness_polys_openings_at_dilations){
            encoding.push(el[2])
        }
        for(let el of this.gate_setup_openings_at_z){
            encoding.push(el[2])
        }
        for(let el of this.gate_selectors_openings_at_z){
            encoding.push(el[1])
        }
        encoding.push(...this.copy_permutation_polys_openings_at_z);
        encoding.push(this.copy_permutation_grand_product_opening_at_z_omega);
        if(this.lookup_s_poly_opening_at_z_omega){
            encoding.push(this.lookup_s_poly_opening_at_z_omega)
        }
        if(this.lookup_grand_product_opening_at_z_omega){
            encoding.push(this.lookup_grand_product_opening_at_z_omega)
        }
        if(this.lookup_t_poly_opening_at_z){
            encoding.push(this.lookup_t_poly_opening_at_z)
        }
        if(this.lookup_t_poly_opening_at_z_omega){
            encoding.push(this.lookup_t_poly_opening_at_z_omega)
        }
        if(this.lookup_selector_poly_opening_at_z){
            encoding.push(this.lookup_selector_poly_opening_at_z)
        }
        if(this.lookup_table_type_poly_opening_at_z){
            encoding.push(this.lookup_table_type_poly_opening_at_z)
        }
        encoding.push(this.quotient_poly_opening_at_z);
        encoding.push(this.linearization_poly_opening_at_z);
        encoding.push(...this.opening_proof_at_z.encode());
        encoding.push(...this.opening_proof_at_z_omega.encode());
        
        return [this.inputs, encoding];
    }
}

class Buf{
    inner: Buffer;

    constructor(buf: Buffer){
        this.inner = buf;
    }

    read_optional_flag(): boolean {
        if(this.read_u64() == BigNumber.from(0)){
            return false
        }else{
            return true
        }
    }

    read_u64(): BigNumber{
        let current = this.inner.slice(0,8);
        this.inner = this.inner.slice(8,this.inner.length);
        return BigNumber.from(current);
    }

    read_fr(): BigNumber {
        let current = this.inner.slice(0,32);
        this.inner = this.inner.slice(32,this.inner.length);
        return BigNumber.from(current);
    }

    read_optional_fr(): BigNumber | undefined {
        if(this.read_optional_flag()){
            return this.read_fr()
        }
    }

    read_curve_affine(): Point {
        const x = BigNumber.from(this.read_fr());
        const y = BigNumber.from(this.read_fr());
        return new Point(x, y);
    } 
    
    read_fr_vec() : BigNumber[]{
        const num_elements = this.read_u64().toNumber();
        let elements : BigNumber[] = [];
        for(let i = 0; i< num_elements; i++){
            elements.push(this.read_fr())
        }
    
        return elements;
    }

    read_curve_affine_vector() : Point[]{
        const num_elements = this.read_u64().toNumber();
        let elements : Point[] = [];
        for(let i = 0; i< num_elements; i++){
            elements.push(this.read_curve_affine())
        }
        return elements;
    }

    read_optional_curve_affine(): Point | undefined{
        if(this.read_optional_flag()){
            return this.read_curve_affine();
        }
    }

    read_tuple_with_one_index(): [BigNumber, BigNumber]{
        const index = this.read_u64();
        const fr = this.read_fr();
        return [index, fr];
    }

    read_tuple_with_one_index_vec(): [BigNumber, BigNumber][]{
        const num_elements = this.read_u64().toNumber();
        let elements = [];
        for(let i = 0; i< num_elements; i++){
            elements.push(this.read_tuple_with_one_index())
        }
        return elements;
    }

    read_tuple_with_two_index(): [BigNumber,  BigNumber, BigNumber]{
        const index0 = this.read_u64();
        const index1 = this.read_u64();
        const fr = this.read_fr();
        return [index0, index1, fr];
    }

    read_tuple_with_two_indexes_vec(): [BigNumber,  BigNumber, BigNumber][]{
        const num_elements = this.read_u64().toNumber();
        let elements = [];
        for(let i = 0; i< num_elements; i++){
            elements.push(this.read_tuple_with_two_index())
        }
        return elements;
    }
}
