SMT models for finding quasidifferential trails in RECTANGLE, KNOT, Speck and Simon.
Using this software requires Python 3 with [pyboolector](https://boolector.github.io/) installed, and Sage for some of the scripts.

## Contents

- `quasidifferential_transition_matrix.sage`: implementation of the 'divide-and-conquer' algorithm to compute quasidifferential transition matrices.
- `arx`: SMT models for Simon-32 and all variants of speck Speck.
- `spn`: SMT models for Rectangle and KNOT. Supporting SageMath scripts to derive conditions on the key from quasidifferential trails with absolute correlation equal to the probability of a characteristic.
