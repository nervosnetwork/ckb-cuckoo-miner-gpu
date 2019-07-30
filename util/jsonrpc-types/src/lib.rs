mod block_template;
mod blockchain;
mod bytes;
mod proposal_short_id;
mod string;

#[derive(Clone, Default, Serialize, Deserialize, PartialEq, Eq, Hash, Debug)]
pub struct BlockNumber(#[serde(with = "string")] pub ckb_core::BlockNumber);

#[derive(Clone, Default, Serialize, Deserialize, PartialEq, Eq, Hash, Debug)]
pub struct Capacity(#[serde(with = "string")] pub ckb_core::Capacity);

#[derive(Clone, Default, Serialize, Deserialize, PartialEq, Eq, Hash, Debug)]
pub struct Cycle(#[serde(with = "string")] pub ckb_core::Cycle);

#[derive(Clone, Default, Serialize, Deserialize, PartialEq, Eq, Hash, Debug)]
pub struct EpochNumber(#[serde(with = "string")] pub ckb_core::EpochNumber);

#[derive(Clone, Default, Serialize, Deserialize, PartialEq, Eq, Hash, Debug)]
pub struct Version(#[serde(with = "string")] pub ckb_core::Version);

#[derive(Clone, Default, Serialize, Deserialize, PartialEq, Eq, Hash, Debug)]
pub struct Timestamp(#[serde(with = "string")] pub u64);

#[derive(Clone, Default, Serialize, Deserialize, PartialEq, Eq, Hash, Debug)]
pub struct Unsigned(#[serde(with = "string")] pub u64);


pub use self::block_template::{
    BlockTemplate, CellbaseTemplate, TransactionTemplate, UncleTemplate,
};
pub use self::blockchain::{
    Block, BlockView, CellInput, CellOutPoint, CellOutput, Header,
    HeaderView, OutPoint, Script, Seal, Transaction, TransactionView, UncleBlock, UncleBlockView, Witness,
};
pub use self::bytes::JsonBytes;
pub use self::proposal_short_id::ProposalShortId;
pub use jsonrpc_core::types::{error, id, params, request, response, version};
pub use serde_derive::{Deserialize, Serialize};
