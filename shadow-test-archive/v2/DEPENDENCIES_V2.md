# V2 dependencies and environment lock

The live capture and outcome validator require:

- a POSIX/Unix-like filesystem and process environment (`fcntl` locks and
  same-filesystem hard links are used for concurrency-safe immutable records);
- Python 3.10 or newer;
- NumPy;
- MDAnalysis.

The development remeasurement and threshold-selection steps additionally
require pandas. The test suite requires pytest. V2 does **not** require
scikit-learn or joblib at runtime.

Policy, enrollment, arm, event, capture, disposition, and validation directories
must remain on local or shared storage that implements atomic rename, hard-link,
and advisory-lock semantics correctly. Do not place the live artifact directory
on an object-store mount that only approximates POSIX behavior.

Do not copy version numbers from this review environment. Run the strict audit,
policy freeze, live monitor, and validator in the actual project environment.
`freeze-policy` records the exact Python, NumPy, and MDAnalysis versions, and
every capture/validation refuses a mismatch. Runtime implementation files are
also SHA-256-bound into the policy.

A container or lockfile may be added before policy freeze. Changing the frozen
software environment or runtime implementation after the first capture creates
a new policy version and a new prospective cohort.
