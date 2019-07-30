use super::{Worker, WorkerMessage};
use crate::pow_message;
use super::CYCLE_LEN;
use byteorder::{ByteOrder, LittleEndian};
use ckb_core::header::Seal;
use ckb_logger::{debug, error};
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;
use numext_fixed_hash::H256;
use std::thread;
use std::time::{Duration, Instant};


extern "C" {
    pub fn c_solve(output: *mut u32, input: *const u8) -> u32;
}

pub struct CuckooSimple {
    start: bool,
    pow_hash: Option<H256>,
    seal_tx: Sender<(H256, Seal)>,
    worker_rx: Receiver<WorkerMessage>,
    seal_candidates_found: u64,
}

impl CuckooSimple {
    pub fn new(
        seal_tx: Sender<(H256, Seal)>,
        worker_rx: Receiver<WorkerMessage>,
    ) -> Self {
        Self {
            start: true,
            pow_hash: None,
            seal_candidates_found: 0,
            seal_tx,
            worker_rx,
        }
    }

    fn poll_worker_message(&mut self) {
        if let Ok(msg) = self.worker_rx.try_recv() {
            match msg {
                WorkerMessage::NewWork(pow_hash) => {
                    self.pow_hash = Some(pow_hash);
                }
                WorkerMessage::Stop => {
                    self.start = false;
                }
                WorkerMessage::Start => {
                    self.start = true;
                }
            }
        }
    }

    #[inline]
    fn solve(&mut self, pow_hash: &H256, nonce: u64) {
        unsafe {
            let mut output = vec![0u32; CYCLE_LEN];
            let ns = c_solve(output.as_mut_ptr(), pow_message(pow_hash, nonce).as_ptr());
            if ns == 1 {
                let mut proof_u8 = vec![0u8; CYCLE_LEN << 2];
                LittleEndian::write_u32_into(&output, &mut proof_u8);
                let seal = Seal::new(nonce, proof_u8.into());
                debug!(
                    "send new found seal, pow_hash {:x}, seal {:?}",
                    pow_hash, seal
                );
                if let Err(err) = self.seal_tx.send((pow_hash.clone(), seal)) {
                    error!("seal_tx send error {:?}", err);
                }
                self.seal_candidates_found += 1;
            }
        }
    }

}


const STATE_UPDATE_DURATION_MILLIS: u128 = 3000;

impl Worker for CuckooSimple {
    fn run<G: FnMut() -> u64>(&mut self, mut rng: G, progress_bar: ProgressBar) {
        let mut state_update_counter = 0usize;
        let mut start = Instant::now();
        loop {
            self.poll_worker_message();
            if self.start {
                if let Some(pow_hash) = self.pow_hash.clone() {
                    self.solve(&pow_hash, rng());
                    state_update_counter += 1;

                    let elapsed = start.elapsed();
                    if elapsed.as_millis() > STATE_UPDATE_DURATION_MILLIS {
                        let elapsed_nanos: f64 = (elapsed.as_secs() * 1_000_000_000
                            + u64::from(elapsed.subsec_nanos()))
                            as f64
                            / 1_000_000_000.0;
                        progress_bar.set_message(&format!(
                            "gps: {:>10.3} / cycles found: {:>10}",
                            state_update_counter as f64 / elapsed_nanos,
                            self.seal_candidates_found,
                        ));
                        progress_bar.inc(1);
                        state_update_counter = 0;
                        start = Instant::now();
                    }
                }
            } else {
                // reset state and sleep
                state_update_counter = 0;
                start = Instant::now();
                thread::sleep(Duration::from_millis(100));
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ckb_pow::{CuckooEngine, PowEngine};
    use crossbeam_channel::unbounded;
    use proptest::prelude::*;

    // fn _cuckoo_solve(pow_hash: &H256, nonce: u64) -> Result<(), TestCaseError> {
    //     let (seal_tx, seal_rx) = unbounded();
    //     let (_worker_tx, worker_rx) = unbounded();
    //     let cuckoo = Cuckoo::new(6, 8);
    //     let mut worker = CuckooSimple::new(cuckoo.clone(), seal_tx, worker_rx);
    //     worker.solve(pow_hash, nonce);
    //     let engine = CuckooEngine { cuckoo };
    //     while let Ok((pow_hash, seal)) = seal_rx.try_recv() {
    //         let (nonce, proof) = seal.destruct();
    //         let message = pow_message(&pow_hash, nonce);
    //         prop_assert!(engine.verify(0, &message, &proof));
    //     }

    //     Ok(())
    // }

    // proptest! {
    //     #[test]
    //     fn cuckoo_solve(h256 in prop::array::uniform32(0u8..), nonce in any::<u64>()) {
    //         _cuckoo_solve(&H256::from_slice(&h256).unwrap(), nonce)?;
    //     }
    // }

    // #[test]
    // fn test_mesg() {
    //     let nonce: u64 = 13999154882437851188;
    //     let msg = [0u8; 32];
    //     let hh: H256 = msg.into();
    //     let result = blake2b_256(pow_message(&hh, nonce).as_ref());
    //     println!("{:?}", &hh[..]);

    // }
}
