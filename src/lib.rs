mod client;
mod config;
mod miner;
mod worker;

use byteorder::{ByteOrder, LittleEndian};
pub use crate::client::Client;
pub use crate::config::{ClientConfig, MinerConfig, WorkerConfig};
pub use crate::miner::Miner;

use ckb_core::block::{Block, BlockBuilder};
use ckb_jsonrpc_types::BlockTemplate;
use std::convert::From;
use ckb_core::difficulty::difficulty_to_target;
use ckb_hash::blake2b_256;
use numext_fixed_hash::H256;
use numext_fixed_uint::U256;

pub struct Work {
    work_id: u64,
    block: Block,
}

impl From<BlockTemplate> for Work {
    fn from(block_template: BlockTemplate) -> Work {
        let work_id = block_template.work_id.clone();
        let block: BlockBuilder = block_template.into();
        let block = block.build();

        Work {
            work_id: work_id.0,
            block,
        }
    }
}

pub fn verify_proof_difficulty(proof: &[u8], difficulty: &U256) -> bool {
    let proof_hash: H256 = blake2b_256(proof).into();
    proof_hash < difficulty_to_target(difficulty)
}

pub fn pow_message(pow_hash: &H256, nonce: u64) -> [u8; 40] {
    let mut message = [0; 40];
    message[8..40].copy_from_slice(&pow_hash[..]);
    LittleEndian::write_u64(&mut message, nonce);
    message
}