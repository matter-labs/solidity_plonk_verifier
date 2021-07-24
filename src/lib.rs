use handlebars::*;
// use serde_json::{Map, Serializer};
use franklin_crypto::bellman::{CurveAffine, Engine, PrimeField, PrimeFieldRepr, pairing::bn256::{Bn256, Fq, Fr}, plonk::{better_better_cs::{cs::Circuit, setup::VerificationKey}, domains::Domain}};
use serde_json::Map;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;

use serde::ser::{Serialize, SerializeStruct, Serializer};

pub struct Commitment<E: Engine> {
    x: E::Fq,
    y: E::Fq,
}
pub struct FieldElement<F: PrimeField>(F);

fn field_element_to_hex<F: PrimeField>(el: &F) -> String {
    let mut buf = [0u8; 32];
    el.into_repr()
        .write_le(&mut buf[..])
        .expect("consume buffer");

    let mut result = String::from("0x");
    result.push_str(&hex::encode(buf));

    result
}


impl<E: Engine> Commitment<E> {
    fn from_affine_point(point: E::G1Affine) -> Self {
        let (x, y) = point.into_xy_unchecked();
        Self { x, y }
    }
}

impl<F: PrimeField> From<F>  for FieldElement<F> {
    fn from(el: F) -> Self {
        Self(el)
    }
}

impl<E: Engine> Serialize for Commitment<E> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // 3 is the number of fields in the struct.
        let mut state = serializer.serialize_struct("Commitment", 2)?;

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
impl<F: PrimeField> Serialize for FieldElement<F> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        // 3 is the number of fields in the struct.
        let mut state = serializer.serialize_struct("FieldElement", 1)?;
        
        let mut init_str = String::from("PairingsBn254.new_fr");
        init_str.push_str("(");
        init_str.push_str(&field_element_to_hex(&self.0));
        init_str.push_str(")");
        state.serialize_field("fr", &init_str)?;
        // state.serialize_field("el", &field_element_to_hex(&self.0))?;
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

pub fn render_verification_key<E: Engine, C: Circuit<E>>(vk: &VerificationKey<E, C>, render_to_path: &str) {

    let mut map = MapWrapper::new();
    let mut handlebars = Handlebars::new();

    // global constants

    // main gate + custom rescue
    map.insert("NUM_GATES", 2); 
     // qL, qR, qM, qO, qConst, qD, qDnext
    map.insert("NUM_SELECTORS", 7);
    // a, b, c, d
    map.insert("STATE_WIDTH", 4); 

    let domain: Domain<E::Fr> = Domain::new_for_size(vk.n as u64).expect("a domain");
    map.insert("domain_size".into(), domain.size);
    map.insert("domain_generator".into(), FieldElement::from(domain.generator));

    let mut setup_commitments = vec![];

    for cmt in vk.gate_setup_commitments.iter() {
        setup_commitments.push(Commitment::<E>::from_affine_point(cmt.clone()))
    }

    map.insert("setup_commitments", setup_commitments);

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
        render_verification_key(&vk, "./contracts/Verifier.sol");
    }
}
