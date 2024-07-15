use std::ops::Rem;

use borsh::{BorshDeserialize, BorshSerialize};
use generic_array::{ArrayLength, GenericArray};
use grid::Grid;
use typenum::{IsNotEqual, U0, U64};

pub const TWO_TO_SIXTEEEN: usize = 1 << 16;

pub trait NonZeroMultipleOf64: ArrayLength + std::fmt::Debug {}
impl<T: ArrayLength + Rem<U64, Output = U0> + IsNotEqual<U0> + std::fmt::Debug> NonZeroMultipleOf64
    for T
{
}

/// Type aliases

pub type Share<Sz> = GenericArray<u8, Sz>;
pub type MerkleRoot = [u8; 32];
pub(crate) type Square<Sz> = Grid<Share<Sz>>;

/// Wrapped Sha256 Merkle proof
pub struct MerkleProof(pub rs_merkle::MerkleProof<rs_merkle::algorithms::Sha256>);

impl From<rs_merkle::MerkleProof<rs_merkle::algorithms::Sha256>> for MerkleProof {
    fn from(value: rs_merkle::MerkleProof<rs_merkle::algorithms::Sha256>) -> Self {
        Self(value)
    }
}

impl std::fmt::Debug for MerkleProof {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Vec::<String>::fmt(&self.0.proof_hashes_hex(), f)
    }
}

impl PartialEq for MerkleProof {
    fn eq(&self, other: &Self) -> bool {
        self.0.proof_hashes() == other.0.proof_hashes()
    }
}

impl BorshSerialize for MerkleProof {
    /// Serialization structure:
    ///
    /// - proof size
    /// - proof as bytes (from library)
    fn serialize<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        let buf = self.0.to_bytes();
        (buf.len() as u64).serialize(writer)?;
        writer.write(&buf).map(drop)
    }
}

impl BorshDeserialize for MerkleProof {
    /// See `serialize` for serialization structure.
    fn deserialize_reader<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        let len = u64::deserialize_reader(reader)?;
        let mut buf = vec![0u8; len as usize];
        reader.read(&mut buf[..])?;
        let proof = rs_merkle::MerkleProof::<rs_merkle::algorithms::Sha256>::from_bytes(&buf)
            .map_err(|err| std::io::Error::new(std::io::ErrorKind::Other, err))?;
        Ok(Self(proof))
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Axis {
    Row,
    Col,
}

/// Creates a new empty share of size `Sz`.
pub fn new_share<Sz: ArrayLength + std::fmt::Debug>() -> Share<Sz> {
    GenericArray::default()
}

impl std::fmt::Display for Axis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Axis::Row => write!(f, "row"),
            Axis::Col => write!(f, "col"),
        }
    }
}

impl Axis {
    pub fn ortho(&self) -> Self {
        match self {
            Axis::Row => Axis::Col,
            Axis::Col => Axis::Row,
        }
    }
}

#[cfg(test)]
mod test {
    use rs_merkle::algorithms::Sha256;
    use rs_merkle::Hasher;
    use rs_merkle::MerkleTree;

    use super::*;

    #[test]
    fn merkle_proof_serialization_roundtrip() {
        let leaves_data = ["a", "b", "c"];
        let leaves = leaves_data.map(|str| Sha256::hash(str.as_bytes()));
        let tree = MerkleTree::<Sha256>::from_leaves(&leaves);
        let proof = MerkleProof(tree.proof(&[0]));
        let mut buf = Vec::<u8>::new();
        proof.serialize(&mut buf).unwrap();
        let desered = MerkleProof::deserialize(&mut &buf[..]).unwrap();
        assert_eq!(proof, desered)
    }
}
