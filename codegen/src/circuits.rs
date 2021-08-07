use franklin_crypto::{
    bellman::{
        plonk::{
            better_better_cs::{
                cs::{
                    Circuit, ConstraintSystem, Gate, GateInternal,
                    LookupTableApplication,
                    PolyIdentifier, Width4MainGateWithDNext,
                },
                gates::selector_optimized_with_d_next::SelectorOptimizedWidth4MainGateWithDNext,
            },
        },
        Engine, Field, PrimeField, SynthesisError,
    },
    plonk::circuit::{
        allocated_num::{AllocatedNum, Num},
        boolean::Boolean,
        custom_rescue_gate::Rescue5CustomGate,
    },
};

use rescue_poseidon::{circuit_generic_hash, CustomGate, HashParams, RescueParams};

#[derive(Clone, Debug)]
pub struct DummyCircuit<E: Engine> {
    input_witness: Option<E::Fr>,
    _m: std::marker::PhantomData<E>,
}
// circuit should have access to at least one lookup table, and at least have one query form it
// circuit should also have at least one sbox(custom gate) in order to pass a test in setup phase
impl<E: Engine> Circuit<E> for DummyCircuit<E> {
    type MainGate = Width4MainGateWithDNext;

    fn synthesize<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<(), SynthesisError> {
        inner_circuit(cs)?;

        Ok(())
    }

    fn declare_used_gates() -> Result<Vec<Box<dyn GateInternal<E>>>, SynthesisError> {
        Ok(vec![
            Self::MainGate::default().into_internal(),
            Rescue5CustomGate::default().into_internal(),
        ])
    }
}

pub struct SelectorOptimizedDummyCircuit;

impl<E: Engine> Circuit<E> for SelectorOptimizedDummyCircuit {
    type MainGate = SelectorOptimizedWidth4MainGateWithDNext;
    fn synthesize<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<(), SynthesisError> {
        inner_circuit(cs)?;
        Ok(())
    }
    fn declare_used_gates() -> Result<Vec<Box<dyn GateInternal<E>>>, SynthesisError> {
        let gates = vec![
            Self::MainGate::default().into_internal(),
            Rescue5CustomGate::default().into_internal(),
        ];
        Ok(gates)
    }
}

fn inner_circuit<E: Engine, CS: ConstraintSystem<E>>(cs: &mut CS) -> Result<(), SynthesisError> {
    for _ in 0..32 {
        let a = Num::alloc(cs, Some(E::Fr::one()))?;
        let b = Num::alloc(cs, Some(E::Fr::zero()))?;
        let flag = Boolean::alloc(cs, Some(true))?;
        let c = Num::conditionally_select(cs, &flag, &a, &b)?;
        let is_equal = Num::equals(cs, &a, &c)?;

        Boolean::enforce_equal(cs, &is_equal, &Boolean::Constant(true))?;
    }

    // add dummy lookup table queries
    let dummy = CS::get_dummy_variable();

    // need to create a table (any)
    let columns = vec![
        PolyIdentifier::VariablesPolynomial(0),
        PolyIdentifier::VariablesPolynomial(1),
        PolyIdentifier::VariablesPolynomial(2),
    ];
    let range_table = LookupTableApplication::new_range_table_of_width_3(2, columns.clone())?;
    let _range_table_name = range_table.functional_name();

    let xor_table = LookupTableApplication::new_xor_table(2, columns.clone())?;
    let _xor_table_name = xor_table.functional_name();

    let and_table = LookupTableApplication::new_and_table(2, columns)?;
    let and_table_name = and_table.functional_name();

    cs.add_table(range_table)?;
    cs.add_table(xor_table)?;
    cs.add_table(and_table)?;

    let binary_x_value = E::Fr::from_str("3").unwrap();
    let binary_y_value = E::Fr::from_str("1").unwrap();

    let t = AllocatedNum::zero(cs);
    let tt = AllocatedNum::one(cs);
    let ttt = t.mul(cs, &tt)?;
    ttt.inputize(cs)?;

    let binary_x = cs.alloc(|| Ok(binary_x_value))?;

    let binary_y = cs.alloc(|| Ok(binary_y_value))?;

    {
        let table = cs.get_table(&and_table_name)?;
        let num_keys_and_values = table.width();

        let and_result_value = table.query(&[binary_x_value, binary_y_value])?[0];

        let binary_z = cs.alloc(|| Ok(and_result_value))?;

        cs.begin_gates_batch_for_step()?;

        let vars = [binary_x, binary_y, binary_z, dummy];
        cs.allocate_variables_without_gate(&vars, &[])?;

        cs.apply_single_lookup_gate(&vars[..num_keys_and_values], table)?;

        cs.end_gates_batch_for_step()?;
    }

    // make single rescue hash to satisfy gate requirements of declaration
    let mut params = RescueParams::default();
    params.use_custom_gate(CustomGate::QuinticWidth4);

    let elem = Num::alloc(cs, Some(E::Fr::from_str("42").unwrap()))?;
    let _ = circuit_generic_hash::<_, _, _, 2, 3, 2>(cs, &[elem, elem], &params, None)?;

    Ok(())
}
