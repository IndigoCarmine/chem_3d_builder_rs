//! Serialization for tests that drive OpenBabel.
//!
//! `openbabel_rs` takes its global lock once per FFI call, which is enough to
//! keep any single call safe but not a *sequence* of them. OpenBabel keeps its
//! force fields as shared singletons, so `generate_3d()` followed by `energy()`
//! is two independent critical sections with another thread's molecule free to
//! pass through the same force field in between.
//!
//! The app never does this: everything that reaches OpenBabel is gated on
//! `MmState::is_running()`, and a minimization owns the molecule outright while
//! it runs. The test harness does, because it runs tests in parallel by default
//! — which showed up as roughly one run in six failing, and two unrelated tests
//! failing together when it did.
//!
//! So every test that drives OpenBabel takes this guard first.

use std::sync::{Mutex, MutexGuard, OnceLock};

/// Hold this for the duration of a test that calls into OpenBabel.
///
/// Poisoning is ignored on purpose: one test failing while holding the guard
/// says nothing about whether OpenBabel is still usable, and turning that into a
/// cascade of confusing failures in every later test would hide the real one.
pub fn ob_guard() -> MutexGuard<'static, ()> {
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(()))
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}
