use handlebars::*;
// use serde_json::{Map, Serializer};
use franklin_crypto::bellman::{CurveAffine, Engine, PrimeField, PrimeFieldRepr, SynthesisError, pairing::bn256::{Bn256, Fq, Fr, Fq2}, plonk::{better_better_cs::{cs::{Circuit, ConstraintSystem, Width4MainGateWithDNext}, setup::VerificationKey}, domains::Domain}};
use serde_json::Map;
use std::io::{Read, Write};

#[macro_use]
extern crate serde_derive;
extern crate serde_json;

use serde::ser::{Serialize, SerializeStruct, Serializer};


const TEMPLATE_FILE_PATH: &str = "./template/verifier.sol";
const VERIFICATION_KEY_FILE_NAME: &str = "VerificationKey.sol";

pub struct DummyCircuitWithRescueCustomGate;

impl<E: Engine> Circuit<E> for DummyCircuitWithRescueCustomGate {
    type MainGate = Width4MainGateWithDNext;
    fn synthesize<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<(), SynthesisError> {
        todo!()
    }
}


pub struct FieldElement(Fr);


impl From<Fr>  for FieldElement {
    fn from(el: Fr) -> Self {
        Self(el)
    }
}

impl Serialize for FieldElement {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("FieldElement", 1)?;
        
        let mut init_str = String::from("PairingsBn254.new_fr");
        init_str.push_str("(");
        init_str.push_str(&field_element_to_hex(&self.0));
        init_str.push_str(")");
        state.serialize_field("el", &init_str)?;        
        state.end()
    }
}

fn field_element_to_hex<F: PrimeField>(el: &F) -> String {
    let mut buf = [0u8; 32];
    el.into_repr()
        .write_be(&mut buf[..])
        .expect("consume buffer");

    let mut result = String::from("0x");
    result.push_str(&hex::encode(buf));

    result
}


pub struct G1Point {
    x: Fq,
    y: Fq,
}

impl G1Point {
    fn from_affine_point(point: <Bn256 as Engine>::G1Affine) -> Self {
        let (x, y) = point.into_xy_unchecked();
        Self { x, y }
    }
}

impl Serialize for G1Point {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("G1Point", 2)?;

        let mut init_str = String::from("PairingsBn254.new_g1");
        init_str.push_str("(");
        init_str.push_str(&field_element_to_hex(&self.x));
        init_str.push_str(",");
        init_str.push_str(&field_element_to_hex(&self.y));
        init_str.push_str(")");
        state.serialize_field("g1", &init_str)?;
        state.end()
    }
}

pub struct G2Point {
    x: Fq2,
    y: Fq2,
}

impl G2Point {
    fn from_affine_point(point: <Bn256 as Engine>::G2Affine) -> Self {
        let (x, y) = point.into_xy_unchecked();
        Self { x, y }
    }
}

impl Serialize for G2Point {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        
        let mut state = serializer.serialize_struct("G2Point", 2)?;
        let mut init_str = String::from("PairingsBn254.new_g2");
        // Encoding of field elements is: X[0] * z + X[1]
        init_str.push_str("([");
        init_str.push_str(&field_element_to_hex(&self.x.c1));
        init_str.push_str(",");
        init_str.push_str(&field_element_to_hex(&self.x.c0));
        init_str.push_str("],[");
        init_str.push_str(&field_element_to_hex(&self.y.c1));
        init_str.push_str(",");
        init_str.push_str(&field_element_to_hex(&self.y.c0));
        init_str.push_str("])");
        state.serialize_field("g2", &init_str)?;
        state.end()
    }
}


struct MapWrapper{
    inner: Map<String, JsonValue>,
}
impl MapWrapper{
    fn new() -> Self{
        Self{
            inner: Map::new()
        }
    }

    fn insert<T: Serialize>(&mut self, key: &str, value: T) -> Option<JsonValue>{
        self.inner.insert(key.into(), to_json(value))
    }
}

/// Read Verification Key from `vk_path`, render solidity verifier 
/// and write its contents into file in `output_path` 
pub fn generate(vk_path: &str, output_path: &str, template_file_path: Option<&str>){  
    let mut reader = std::io::BufReader::new(
        std::fs::File::open(vk_path).expect("should open file"),
    );
    let mut final_output_path = String::new();
    final_output_path.push_str(output_path);
    final_output_path.push_str("/");
    final_output_path.push_str(VERIFICATION_KEY_FILE_NAME);
    println!("final output {}", final_output_path);

    let mut writer = std::io::BufWriter::new(
        std::fs::File::create(final_output_path).expect("should open output file")
    );

    let template_file_path = if let Some(path) = template_file_path{
        path
    }else{
        TEMPLATE_FILE_PATH
    };


    render_verification_key::<DummyCircuitWithRescueCustomGate, _, _>(&template_file_path, &mut reader, &mut writer)
}

fn render_verification_key<C: Circuit<Bn256>, R: Read, W: Write>(template_file_path: &str, reader: &mut R, writer: &mut W) {
    let vk = VerificationKey::<Bn256, DummyCircuitWithRescueCustomGate>::read(reader).expect("parsed vk");

    let mut map = MapWrapper::new();
    let mut handlebars = Handlebars::new();

    // global constants

    // main gate + custom rescue
    map.insert("NUM_GATES", vk.gate_selectors_commitments.len()); 
     // qL, qR, qM, qO, qConst, qD, qDnext
    map.insert("NUM_SELECTORS", vk.gate_setup_commitments.len());
    // a, b, c, d
    map.insert("STATE_WIDTH", vk.state_width); 
    map.insert("DNEXT_INDEX", vk.state_width-1); 
    map.insert("NUM_G2_ELS", vk.g2_elements.len());     
    map.insert("NUM_LOOKUP_TABLES", vk.lookup_tables_commitments.len()); 
    map.insert("SERIALIZED_PROOF_LENGTH", 44);  // TODO calculat length
    map.insert("NUM_ALPHA_CHALLENGES", 8);  // TODO

    // domain
    map.insert("num_inputs".into(), vk.num_inputs);
    let domain: Domain<Fr> = Domain::new_for_size(vk.n as u64).expect("a domain");
    map.insert("domain_size".into(), domain.size);
    map.insert("domain_generator".into(), FieldElement::from(domain.generator));

    // G1Points
    let mut gate_setup_commitments = vec![];
    for cmt in vk.gate_setup_commitments.iter() {
        gate_setup_commitments.push(G1Point::from_affine_point(cmt.clone()))
    }
    map.insert("gate_setup_commitments", gate_setup_commitments);

    let mut gate_selectors_commitments = vec![];
    for cmt in vk.gate_selectors_commitments.iter() {
        gate_selectors_commitments.push(G1Point::from_affine_point(cmt.clone()))
    }
    map.insert("gate_selectors_commitments", gate_selectors_commitments);
    
    let mut permutation_commitments = vec![];
    for cmt in vk.permutation_commitments.iter() {
        permutation_commitments.push(G1Point::from_affine_point(cmt.clone()))
    }
    map.insert("permutation_commitments", permutation_commitments);    
    
    map.insert("total_lookup_entries_length", vk.total_lookup_entries_length);    
    map.insert("lookup_selector_commitment", G1Point::from_affine_point(vk.lookup_selector_commitment.unwrap_or(<Bn256 as Engine>::G1Affine::zero())));
    if vk.total_lookup_entries_length > 0 {
        assert!(vk.lookup_selector_commitment.is_some());
    }
    let mut lookup_tables_commitments = vec![];
    for cmt in vk.lookup_tables_commitments.iter() {
        lookup_tables_commitments.push(G1Point::from_affine_point(cmt.clone()))
    }
    map.insert("lookup_tables_commitments", lookup_tables_commitments);
    map.insert("lookup_table_type_commitment", G1Point::from_affine_point(vk.lookup_table_type_commitment.unwrap_or(<Bn256 as Engine>::G1Affine::zero())));
    
    // non residues
    let mut non_residues = vec![];
    for el in vk.non_residues.iter() {
        non_residues.push(FieldElement::from(el.clone()));
    }
    map.insert("non_residues", non_residues);
    
    // pairing g2 elements
    let mut g2_elements = vec![];
    for point in vk.g2_elements.iter() {
        g2_elements.push(G2Point::from_affine_point(point.clone()));
    }
    map.insert("g2_elements", g2_elements);

    // register template from a file and assign a name to it
    handlebars
        .register_template_file("contract", template_file_path)
        .expect("must read the template");


    let rendered = handlebars.render("contract", &map.inner).unwrap();
    
    writer
        .write(rendered.as_bytes())
        .expect("must write to file");
}

#[cfg(test)]
mod tests {
    use super::*;
    use franklin_crypto::bellman::{Engine, SynthesisError, pairing::bn256::*, plonk::better_better_cs::{cs::{Circuit, ConstraintSystem, Width4MainGateWithDNext}}};

    #[test]
    fn test_render_verification_key() {
        generate("./../block_vk_20_keccak.key", "./../hardhat/contracts", None);
    }
}
