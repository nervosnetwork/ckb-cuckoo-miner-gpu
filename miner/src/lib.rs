mod client;
mod config;
mod miner;
mod worker;

pub use crate::client::Client;
pub use crate::config::MinerConfig;
pub use crate::miner::Miner;

use ckb_core::block::{Block, BlockBuilder};
use ckb_jsonrpc_types::BlockTemplate;
use std::convert::From;

#[derive(Debug, Clone)]
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
