use franklin_crypto::bellman::plonk::better_better_cs::gates::selector_optimized_with_d_next::SelectorOptimizedWidth4MainGateWithDNext;
use franklin_crypto::bellman::{
    pairing::bn256::{Bn256, Fr},
    plonk::{
        better_better_cs::{
            cs::{Circuit, Width4MainGateWithDNext},
            setup::VerificationKey,
        },
        domains::Domain,
    },
    CurveAffine, Engine,
};
use handlebars::*;
use serde_json::Map;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::str::FromStr;

use serde::ser::Serialize;

use crate::{
    circuits::{DummyCircuit, SelectorOptimizedDummyCircuit},
    serialize::{FieldElement, G1Point, G2Point},
};

const TEMPLATE_FILE_PATH: &str = "./template/verifier.sol";
const VERIFICATION_KEY_FILE_NAME: &str = "VerificationKey.sol";

#[derive(Clone, Debug)]
pub enum MainGateType {
    Standard,
    SelectorOptimized,
}

#[derive(Clone, Debug)]
pub struct UnknownMainGate;

impl ToString for UnknownMainGate {
    fn to_string(&self) -> String {
        String::from("Unkown main gate! ")
    }
}

impl FromStr for MainGateType {
    type Err = UnknownMainGate;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s == "std_width4" {
            return Ok(Self::Standard);
        } else if s == "selector_optimized_width4" {
            return Ok(Self::SelectorOptimized);
        } else {
            Err(UnknownMainGate)
        }
    }
}
/// Reads a VerificationKey from a file that belongs to a circuit
/// which declares a standard width4 or selector optimized main gate
/// and writes rendered solidity file into a file in `output_path`
pub fn generate(
    main_gate_type: MainGateType,
    vk_path: PathBuf,
    output_path: PathBuf,
    template_file_path: Option<&str>,
) {
    match main_gate_type {
        MainGateType::Standard => {
            inner_generate::<DummyCircuit<Bn256>>(vk_path, output_path, template_file_path)
        }
        MainGateType::SelectorOptimized => inner_generate::<SelectorOptimizedDummyCircuit>(
            vk_path,
            output_path,
            template_file_path,
        ),
    }
}

struct MapWrapper {
    inner: Map<String, JsonValue>,
}
impl MapWrapper {
    fn new() -> Self {
        Self { inner: Map::new() }
    }

    fn insert<T: Serialize>(&mut self, key: &str, value: T) -> Option<JsonValue> {
        self.inner.insert(key.into(), to_json(value))
    }
}

const EXISTING_GATES: [GateType; 3] = [
    GateType::MainGate("main gate of width 4 with D_next"),
    GateType::MainGate("main gate of width 4 with D_next and selector optimization"),
    GateType::CustomGate("Alpha=5 custom gate for Rescue/Poseidon"),
];

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum GateType<'a> {
    MainGate(&'a str),
    CustomGate(&'a str),
}

impl<'a> GateType<'a> {
    fn from_name(name: &str) -> Option<Self> {
        let mut found_gate = None;
        for gate in EXISTING_GATES.iter() {
            match gate {
                GateType::MainGate(gate_name) | GateType::CustomGate(gate_name) => {
                    if *gate_name == name {
                        found_gate = Some(gate.clone());
                    }
                }
            }
        }
        found_gate
    }

    fn is_main_gate(&self) -> bool {
        match self {
            Self::MainGate(_) => true,
            _ => false,
        }
    }
}

fn inner_generate<C: Circuit<Bn256>>(
    vk_path: PathBuf,
    mut output_path: PathBuf,
    template_file_path: Option<&str>,
) {
    let mut reader =
        std::io::BufReader::new(std::fs::File::open(vk_path).expect("should open file"));
    output_path.push(VERIFICATION_KEY_FILE_NAME);

    let mut writer = std::io::BufWriter::new(
        std::fs::File::create(output_path).expect("should open output file"),
    );

    let template_file_path = if let Some(path) = template_file_path {
        path
    } else {
        TEMPLATE_FILE_PATH
    };

    render_verification_key::<C, _, _>(&template_file_path, &mut reader, &mut writer)
}

fn render_verification_key<C: Circuit<Bn256>, R: Read, W: Write>(
    template_file_path: &str,
    reader: &mut R,
    writer: &mut W,
) {
    let vk = VerificationKey::<Bn256, C>::read(reader).expect("parsed vk");

    let mut map = MapWrapper::new();
    let mut handlebars = Handlebars::new();

    let declared_gates = C::declare_used_gates().unwrap();
    let num_gates = EXISTING_GATES.len() - 1;
    assert!(declared_gates.len() == num_gates);
    // TODO
    let mut main_gate = None;
    for gate in declared_gates {
        let new_gate = GateType::from_name(gate.name()).expect("a already defined gate");

        if new_gate.is_main_gate() {
            main_gate = Some(new_gate)
        }
    }

    let main_gate = main_gate.expect("circuit should contain a main gate");

    let (num_main_gate_selectors, ab_coeff_idx, constant_coeff_idx, d_next_coeff_idx, ac_coeff_idx) =
        if main_gate == EXISTING_GATES[0] {
            (
                7,
                Width4MainGateWithDNext::AB_MULTIPLICATION_TERM_COEFF_INDEX,
                Width4MainGateWithDNext::CONSTANT_TERM_COEFF_INDEX,
                Width4MainGateWithDNext::D_NEXT_TERM_COEFF_INDEX,
                None,
            )
        } else if main_gate == EXISTING_GATES[1] {
            (
                8,
                SelectorOptimizedWidth4MainGateWithDNext::AB_MULTIPLICATION_TERM_COEFF_INDEX,
                SelectorOptimizedWidth4MainGateWithDNext::CONSTANT_TERM_COEFF_INDEX,
                SelectorOptimizedWidth4MainGateWithDNext::D_NEXT_TERM_COEFF_INDEX,
                Some(SelectorOptimizedWidth4MainGateWithDNext::AC_MULTIPLICATION_TERM_COEFF_INDEX),
            )
        } else {
            unimplemented!("unknown gate");
        };

    let is_selector_optimized_main_gate = ac_coeff_idx.is_some();
    map.insert(
        "is_selector_optimized_main_gate",
        is_selector_optimized_main_gate,
    );

    // main gate + custom rescue
    assert_eq!(vk.gate_selectors_commitments.len(), num_gates);
    map.insert("NUM_GATES", vk.gate_selectors_commitments.len());
    map.insert("MAIN_GATE_AB_COEFF_IDX", ab_coeff_idx);
    map.insert("CONSTANT_TERM_COEFF_INDEX", constant_coeff_idx);
    map.insert("D_NEXT_TERM_COEFF_INDEX", d_next_coeff_idx);
    assert_eq!(ab_coeff_idx, 4);
    if is_selector_optimized_main_gate {
        assert!(ac_coeff_idx.unwrap() != 0, "AC coeff can't be zero");
        map.insert(
            "MAIN_GATE_AC_COEFF_IDX",
            ac_coeff_idx.expect("AC coeff idx if selector optimized main gate"),
        );
    }
    assert_eq!(vk.gate_setup_commitments.len(), num_main_gate_selectors);
    map.insert("NUM_MAIN_GATE_SELECTORS", num_main_gate_selectors);
    // a, b, c, d
    map.insert("STATE_WIDTH", vk.state_width);
    map.insert("DNEXT_INDEX", vk.state_width - 1);
    map.insert("NUM_G2_ELS", vk.g2_elements.len());
    map.insert("NUM_LOOKUP_TABLES", vk.lookup_tables_commitments.len());
    map.insert("SERIALIZED_PROOF_LENGTH", 44); // TODO calculate length
    map.insert("NUM_ALPHA_CHALLENGES", 8); // TODO

    // domain
    map.insert("num_inputs".into(), vk.num_inputs);
    let domain: Domain<Fr> = Domain::new_for_size(vk.n as u64).expect("a domain");
    map.insert("domain_size".into(), domain.size);
    map.insert(
        "domain_generator".into(),
        FieldElement::from(domain.generator),
    );

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

    map.insert(
        "total_lookup_entries_length",
        vk.total_lookup_entries_length,
    );
    map.insert(
        "lookup_selector_commitment",
        G1Point::from_affine_point(
            vk.lookup_selector_commitment
                .unwrap_or(<Bn256 as Engine>::G1Affine::zero()),
        ),
    );
    if vk.total_lookup_entries_length > 0 {
        assert!(vk.lookup_selector_commitment.is_some());
    }
    let mut lookup_tables_commitments = vec![];
    for cmt in vk.lookup_tables_commitments.iter() {
        lookup_tables_commitments.push(G1Point::from_affine_point(cmt.clone()))
    }
    map.insert("lookup_tables_commitments", lookup_tables_commitments);
    map.insert(
        "lookup_table_type_commitment",
        G1Point::from_affine_point(
            vk.lookup_table_type_commitment
                .unwrap_or(<Bn256 as Engine>::G1Affine::zero()),
        ),
    );

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
        .expect(&format!(
            "must read the template at path {}",
            template_file_path
        ));

    let rendered = handlebars.render("contract", &map.inner).unwrap();

    writer
        .write(rendered.as_bytes())
        .expect("must write to file");
}
