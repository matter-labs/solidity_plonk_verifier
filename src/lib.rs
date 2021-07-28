use handlebars::*;
// use serde_json::{Map, Serializer};
use franklin_crypto::bellman::{CurveAffine, Engine, PrimeField, PrimeFieldRepr, pairing::bn256::{Bn256, Fq, Fr, Fq2}, plonk::{better_better_cs::{cs::Circuit, setup::VerificationKey}, domains::Domain}};
use serde_json::Map;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;

use serde::ser::{Serialize, SerializeStruct, Serializer};



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
        .write_le(&mut buf[..])
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
        // 3 is the number of fields in the struct.
        let mut state = serializer.serialize_struct("G2Point", 2)?;
        let mut init_str = String::from("PairingsBn254.new_g2");

        init_str.push_str("([");
        init_str.push_str(&field_element_to_hex(&self.x.c0));
        init_str.push_str(",");
        init_str.push_str(&field_element_to_hex(&self.x.c1));
        init_str.push_str("],[");
        init_str.push_str(&field_element_to_hex(&self.y.c0));
        init_str.push_str(",");
        init_str.push_str(&field_element_to_hex(&self.y.c1));
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

pub fn render_verification_key<C: Circuit<Bn256>>(vk: &VerificationKey<Bn256, C>, render_to_path: &str) {

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
    map.insert("SERIALIZED_PROOF_LENGTH", 0); 

    // domain
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
        .register_template_file("contract", "./template/verifier.sol")
        .expect("must read the template");

    let mut writer =
        std::io::BufWriter::with_capacity(1 << 24, std::fs::File::create(render_to_path).unwrap());

    let rendered = handlebars.render("contract", &map.inner).unwrap();

    use std::io::Write;
    writer
        .write(rendered.as_bytes())
        .expect("must write to file");
}

#[cfg(test)]
mod tests {
    use super::*;
    use franklin_crypto::bellman::{Engine, SynthesisError, pairing::bn256::*, plonk::better_better_cs::{cs::{Circuit, ConstraintSystem, Width4MainGateWithDNext}}};

    pub struct TestCircuit;

    impl<E: Engine> Circuit<E> for TestCircuit {
        type MainGate = Width4MainGateWithDNext;
        fn synthesize<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<(), SynthesisError> {
            todo!()
        }
    }

    #[test]
    fn test_render_verification_key() {
        let mut reader = std::io::BufReader::with_capacity(
            1 << 24,
            // std::fs::File::open("./xor_vk.key").unwrap()
            std::fs::File::open("./vk_17.key").unwrap(),
        );

        let vk = VerificationKey::<Bn256, TestCircuit>::read(&mut reader).expect("parsed vk");
        render_verification_key(&vk, "./hardhat/contracts/VerificationKey.sol");
    }
}
