# Arkworks Adapter for Lambdaworks Groth16 Backend

This crate enables circuits written in Arkworks to be proven / verified via Lambdaworks Groth16 backend.

Crate exposes [to_lambda](./src/lib.rs#to_lambda) function for <b><span style="color: #57a14d">ConstraintSystemRef</span></b> type of Arkworks, which is expected to carry all constraints, witnesses, and public variable assignments:

```rust
pub fn to_lambda<F: PrimeField>(cs: &ConstraintSystemRef<F>) -> (QuadraticArithmeticProgram, Vec<FrElement>)
```

It returns a Lambdaworks-compatible QAP struct alongside variable assignments. Please note that public variable assignments are bundled with witnesses, and this vector of field elements is called **witness** alltogether.

```rust
let (qap, w) = to_lambda(&cs);
```

After this point, typical steps of Groth16 can be performed using Lamdaworks: setup, prove, verify

```rust
let (pk, vk) = setup(&qap);

let proof = Prover::prove(&w, &qap, &pk);

let public_inputs = &w[..qap.num_of_public_inputs];
let accept = verify(&vk, &proof, public_inputs);

assert!(accept);
```

## Full Example

A linear exponentiation example on BLS12-381 can be found here.
Please check [integration_tests.rs](./src/integration_tests.rs) for more examples.

```rust

use crate::to_lambda;
use ark_bls12_381::Fr;
use ark_relations::{lc, r1cs::ConstraintSystem, r1cs::Variable};
use lambdaworks_groth16::{setup, verify, Prover};
use rand::Rng;

// ...
// ...

let mut rng = rand::thread_rng();
let x = rng.gen::<u64>();
let exp = rng.gen::<u8>();

// Define the circuit using Arkworks

let cs = ConstraintSystem::<Fr>::new_ref();

let x = Fr::from(x);
let mut _x = cs.new_witness_variable(|| Ok(x)).unwrap();

let mut acc = Fr::from(x);
let mut _acc = cs.new_witness_variable(|| Ok(x)).unwrap();

for _ in 0..exp - 1 {
	acc *= x;
	let _new_acc = cs.new_witness_variable(|| Ok(acc)).unwrap();
	cs.enforce_constraint(lc!() + _acc, lc!() + _x, lc!() + _new_acc)
		.unwrap();
	_acc = _new_acc;
}

let _out = cs.new_input_variable(|| Ok(acc)).unwrap();
cs.enforce_constraint(lc!() + _out, lc!() + Variable::One, lc!() + _acc)
	.unwrap();

// Make Lambdaworks-compatible
let (qap, w) = to_lambda(&cs);

// Use Lambdaworks Groth16 backend

let (pk, vk) = setup(&qap);

let proof = Prover::prove(&w, &qap, &pk);

let public_inputs = &w[..qap.num_of_public_inputs];
let accept = verify(&vk, &proof, public_inputs);
assert!(accept);
```
