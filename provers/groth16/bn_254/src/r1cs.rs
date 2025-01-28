use crate::common::FrElement;
use lambdaworks_math::field::{element::FieldElement, traits::IsField};

// To be improved with a front-end implementation
// TODO: Use CS in Groth16 tests instead of a plain QAP
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConstraintSystem<F: IsField> {
    pub constraints: R1CS,
    pub witness: Vec<FieldElement<F>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Constraint {
    pub a: Vec<FrElement>,
    pub b: Vec<FrElement>,
    pub c: Vec<FrElement>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct R1CS {
    pub constraints: Vec<Constraint>,
    pub number_of_inputs: usize,
}

impl R1CS {
    pub fn from_matrices(
        a: Vec<Vec<FrElement>>,
        b: Vec<Vec<FrElement>>,
        c: Vec<Vec<FrElement>>,
        number_of_inputs: usize,
    ) -> Self {
        Self {
            constraints: (0..a.len())
                .map(|i| Constraint {
                    a: a[i].clone(),
                    b: b[i].clone(),
                    c: c[i].clone(),
                })
                .collect(),
            number_of_inputs,
        }
    }

    pub fn number_of_constraints(&self) -> usize {
        self.constraints.len()
    }

    pub fn witness_size(&self) -> usize {
        self.constraints[0].a.len()
    }
}
