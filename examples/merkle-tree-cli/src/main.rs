mod commands;
use clap::Parser;
use commands::{MerkleArgs, MerkleEntity};
use lambdaworks_crypto::hash::poseidon::starknet::PoseidonCairoStark252;
use lambdaworks_crypto::merkle_tree::{
    backends::field_element::TreePoseidon, merkle::MerkleTree, proof::Proof,
};
use lambdaworks_math::field::{
    element::FieldElement, fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
};
use std::io::BufWriter;
use std::{
    fs::{self, File},
    io::{self, Write},
};

type FE = FieldElement<Stark252PrimeField>;

fn load_fe_from_file(file_path: &String) -> Result<FE, io::Error> {
    FE::from_hex(&fs::read_to_string(file_path)?.replace('\n', ""))
        .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("{:?}", e)))
}

fn load_tree_values(tree_path: &String) -> Result<Vec<FE>, io::Error> {
    Ok(fs::read_to_string(tree_path)?
        .split(';')
        .map(FE::from_hex_unchecked)
        .collect())
}

fn generate_merkle_tree(tree_path: String) -> Result<(), io::Error> {
    let values: Vec<FE> = load_tree_values(&tree_path)?;

    let merkle_tree = MerkleTree::<TreePoseidon<PoseidonCairoStark252>>::build(&values)
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "requested empty tree"))?;
    let root = merkle_tree.root.representative().to_string();
    println!("Generated merkle tree with root: {:?}", root);

    let generated_tree_path = tree_path.replace(".csv", ".json");
    let file = File::create(generated_tree_path)?;
    let mut writer = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, &merkle_tree)?;
    println!("Saved tree to file");
    Ok(())
}

fn generate_merkle_proof(tree_path: String, pos: usize) -> Result<(), io::Error> {
    let values: Vec<FE> = load_tree_values(&tree_path)?;
    let merkle_tree = MerkleTree::<TreePoseidon<PoseidonCairoStark252>>::build(&values)
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "requested empty tree"))?;

    let Some(proof) = merkle_tree.get_proof(pos) else {
        return Err(io::Error::new(io::ErrorKind::Other, "Index out of bounds"));
    };

    let proof_path = tree_path.replace(".csv", format!("_proof_{pos}.json").as_str());
    let file = File::create(proof_path)?;
    let mut writer = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, &proof)?;
    writer.flush()
}

fn verify_merkle_proof(
    root_path: String,
    index: usize,
    proof_path: String,
    leaf_path: String,
) -> Result<(), io::Error> {
    let root_hash: FE = load_fe_from_file(&root_path)?;

    let leaf: FE = load_fe_from_file(&leaf_path)?;

    let file_str = fs::read_to_string(proof_path)?;
    let proof: Proof<FE> = serde_json::from_str(&file_str)?;

    match proof.verify::<TreePoseidon<PoseidonCairoStark252>>(&root_hash, index, &leaf) {
        true => println!("\x1b[32mMerkle proof verified succesfully\x1b[0m"),
        false => println!("\x1b[31mMerkle proof failed verifying\x1b[0m"),
    }

    Ok(())
}

fn main() {
    let args: MerkleArgs = MerkleArgs::parse();
    if let Err(e) = match args.entity {
        MerkleEntity::GenerateTree(args) => generate_merkle_tree(args.tree_path),
        MerkleEntity::GenerateProof(args) => generate_merkle_proof(args.tree_path, args.position),
        MerkleEntity::VerifyProof(args) => {
            verify_merkle_proof(args.root_path, args.index, args.proof_path, args.leaf_path)
        }
    } {
        println!("Error while running command: {:?}", e);
    }
}
