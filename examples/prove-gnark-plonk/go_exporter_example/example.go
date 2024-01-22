// Copyright 2020 ConsenSys AG
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package main

import (
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"reflect"
	"strconv"
	"unsafe"

	"encoding/json"

	"github.com/consensys/gnark-crypto/ecc"
	fr_bls12381 "github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
	"github.com/consensys/gnark/backend/plonk"
	plonk_bls12381 "github.com/consensys/gnark/backend/plonk/bls12-381"
	"github.com/consensys/gnark/backend/witness"
	cs "github.com/consensys/gnark/constraint/bls12-381"
	"github.com/consensys/gnark/frontend/cs/scs"
	"github.com/consensys/gnark/test"

	"github.com/consensys/gnark/frontend"
)

type SerializedCircuit struct {
	N           int
	N_Padded    uint64
	Omega       string
	Input       []string
	Ql          []string
	Qr          []string
	Qm          []string
	Qo          []string
	Qc          []string
	A           []string
	B           []string
	C           []string
	Permutation []int64
}

func witness_to_json_array(fullWitness witness.Witness) string {

	witnessVec := fullWitness.Vector()
	witnessVecString := fmt.Sprint(witnessVec)
	inputJSON := fmt.Sprint("{\"witness\": ", witnessVecString, "}\n")

	fmt.Println("Witness json")
	fmt.Println(inputJSON)

	// Unmarshal JSON into a map
	var data map[string]interface{}
	json.Unmarshal([]byte(inputJSON), &data)

	// Convert array values to strings
	witnessArray, _ := data["witness"].([]interface{})

	var witnessStrings []string
	for _, value := range witnessArray {
		// Convert the value to a string
		strValue := strconv.Itoa(int(value.(float64)))
		witnessStrings = append(witnessStrings, strValue)
	}

	// Update the map with the converted array
	data["witness"] = witnessStrings

	// Marshal the updated map back to JSON
	outputJSON, _ := json.Marshal(data)

	return string(outputJSON)
}

func ToJSON(_r1cs *cs.SparseR1CS, pk *plonk_bls12381.ProvingKey, fullWitness witness.Witness, witnessPublic fr_bls12381.Vector) { // n
	nbConstraints := _r1cs.GetNbConstraints()
	nbPublic := len(_r1cs.Public)
	n := nbConstraints + nbPublic
	omega := pk.Domain[0].Generator.Text(16)

	// Ql, Qm, Qr, Qo, Qk, S1, S2, S3
	var Ql, Qr, Qm, Qo, Qc []string

	for i := 0; i < nbPublic; i++ {
		var minus_one fr_bls12381.Element
		minus_one = fr_bls12381.NewElement(1)
		minus_one.Neg(&minus_one)
		zero := fr_bls12381.NewElement(0)
		Ql = append(Ql, minus_one.Text(16))
		Qr = append(Qr, zero.Text(16))
		Qm = append(Qm, zero.Text(16))
		Qo = append(Qo, zero.Text(16))
		Qc = append(Qc, zero.Text(16))
	}

	constraint_list := _r1cs.GetSparseR1Cs()

	fmt.Println("Coeffs")
	fmt.Println(_r1cs.Coefficients[0].Text(16))
	fmt.Println(_r1cs.Coefficients[1].Text(16))
	fmt.Println(_r1cs.Coefficients[2].Text(16))
	fmt.Println(_r1cs.Coefficients[3].Text(16))
	fmt.Println("")

	for i := 0; i < nbConstraints; i++ { // constraints
		Ql = append(Ql, _r1cs.Coefficients[int(constraint_list[i].QL)].Text(16))

		fmt.Println("Constraint")
		fmt.Println("Index QL,QR,QO, QM, QC")
		fmt.Println(constraint_list[i].QL)
		fmt.Println(constraint_list[i].QR)
		fmt.Println(constraint_list[i].QO)
		fmt.Println(constraint_list[i].QM)
		fmt.Println(constraint_list[i].QC)

		fmt.Println("Values")
		fmt.Println(_r1cs.Coefficients[int(constraint_list[i].QL)].Text(16))
		fmt.Println(_r1cs.Coefficients[int(constraint_list[i].QR)].Text(16))
		fmt.Println(_r1cs.Coefficients[int(constraint_list[i].QO)].Text(16))
		fmt.Println("")

		Qr = append(Qr, _r1cs.Coefficients[int(constraint_list[i].QR)].Text(16))

		var new_Qm fr_bls12381.Element
		/*
			var new_Qm fr_bls12381.Element
			new_Qm.Set(&_r1cs.Coefficients[_r1cs.Constraints[i].M[0].CoeffID()]).Mul(&new_Qm, &_r1cs.Coefficients[_r1cs.Constraints[i].M[1].CoeffID()])

		*/
		new_Qm.Set(&_r1cs.Coefficients[int(constraint_list[i].QM)])
		// new_Qm.Double(&new_Qm)
		Qm = append(Qm, new_Qm.Text(16))

		Qo = append(Qo, _r1cs.Coefficients[int(constraint_list[i].QO)].Text(16))
		Qc = append(Qc, _r1cs.Coefficients[int(constraint_list[i].QC)].Text(16))
	}

	// Witness
	// opt, _ := backend.NewProverConfig()
	var _solution, _ = _r1cs.Solve(fullWitness)
	abc := _solution.(*cs.SparseR1CSSolution)

	fmt.Println("L = A")
	fmt.Println(abc.L)
	fmt.Println("R = B")
	fmt.Println(abc.R)
	fmt.Println("O = C")
	fmt.Println(abc.O)

	var a, b, c []string

	for i := 0; i < abc.L.Len(); i++ {
		a = append(a, abc.L[i].Text(16))
		b = append(b, abc.R[i].Text(16))
		c = append(c, abc.O[i].Text(16))
	}

	var input []string
	for i := 0; i < len(_r1cs.Public); i++ {
		input = append(input, witnessPublic[i].Text(16))
	}

	/*
		Permutation is a private field, and for a reason, they are changing the API a lot here. It's a bit different in the current main.

		This code needs to be update. Currently working with Gnark 9.1
	*/

	rs := reflect.ValueOf(pk).Elem()
	rf := rs.Field(0)
	rf = reflect.NewAt(rf.Type(), unsafe.Pointer(rf.UnsafeAddr())).Elem()

	var trace_pointer plonk_bls12381.Trace = rf.Interface().(plonk_bls12381.Trace)

	s := trace_pointer.S
	fmt.Println("Trace Permutation")
	fmt.Println(s)

	data := SerializedCircuit{
		N:           n,
		Omega:       omega,
		N_Padded:    pk.Domain[0].Cardinality,
		Input:       input,
		Ql:          Ql,
		Qr:          Qr,
		Qm:          Qm,
		Qo:          Qo,
		Qc:          Qc,
		A:           a,
		B:           b,
		C:           c,
		Permutation: s,
	}
	file, _ := json.MarshalIndent(data, "", " ")
	_ = ioutil.WriteFile("frontend_precomputed_values.json", file, 0644)

	witnessVec := fullWitness.Vector()
	println()
	witnessVecString := fmt.Sprint(witnessVec)

	fmt.Println(witnessVecString)

	f, _ := os.Create("witness.json")
	fmt.Fprintf(f, witness_to_json_array(fullWitness))

	// witnessFormatted := fmt.Sprintln(witnessVec)

	// f, _ := os.Create("witness.json")
	// fmt.Fprintln(f, witnessVec)

	// schema, _ := frontend.NewSchema(&fullWitness)
	// json_wit, _ := fullWitness.ToJSON(schema)
	// _ = ioutil.WriteFile("witness.json", witnessFormatted, 0644)

}

// In this example we show how to use PLONK with KZG commitments. The circuit that is
// showed here is the same as in ../exponentiate.

// Circuit y == x**e
// only the bitSize least significant bits of e are used
type Circuit struct {
	// tagging a variable is optional
	// default uses variable name and secret visibility.
	X frontend.Variable `gnark:",public"`
	Y frontend.Variable `gnark:",public"`
}

// Define declares the circuit's constraints
// y == x*x*x
func (circuit *Circuit) Define(api frontend.API) error {
	api.AssertIsEqual(circuit.Y, api.Mul(circuit.X, circuit.X, circuit.X))
	return nil
}

func main() {

	var circuit Circuit

	// building the circuit...
	ccs, err := frontend.Compile(ecc.BLS12_381.ScalarField(), scs.NewBuilder, &circuit)

	if err != nil {
		fmt.Println("circuit compilation error")
	}

	// create the necessary data for KZG.
	// This is a toy example, normally the trusted setup to build ZKG
	// has been ran before.
	// The size of the data in KZG should be the closest power of 2 bounding //
	// above max(nbConstraints, nbVariables).
	_r1cs := ccs.(*cs.SparseR1CS)
	srs, err := test.NewKZGSRS(_r1cs)
	if err != nil {
		panic(err)
	}

	// Correct data: the proof passes
	{
		// Witnesses instantiation. Witness is known only by the prover,
		// while public w is a public data known by the verifier.
		var w Circuit
		w.X = 3
		w.Y = 27

		witnessFull, err := frontend.NewWitness(&w, ecc.BLS12_381.ScalarField())

		if err != nil {
			log.Fatal(err)
		}

		witnessPublic, err := frontend.NewWitness(&w, ecc.BLS12_381.ScalarField(), frontend.PublicOnly())
		if err != nil {
			log.Fatal(err)
		}

		// public data consists the polynomials describing the constants involved
		// in the constraints, the polynomial describing the permutation ("grand
		// product argument"), and the FFT domains.
		pk, _, err := plonk.Setup(ccs, srs)

		// cs implements io.WriterTo
		// var buf bytes.Buffer
		// ccs.WriteTo(&buf)

		// a, err := json.Marshal(pk)
		// fmt.Println(buf)
		// fmt.Println(string(buf.Bytes()))

		//_, err := plonk.Setup(r1cs, kate, &publicWitness)
		if err != nil {
			log.Fatal(err)
		}

		publicWitness, _ := witnessPublic.Vector().(fr_bls12381.Vector)

		ToJSON(_r1cs, pk.(*plonk_bls12381.ProvingKey), witnessFull, publicWitness)
	}
}
