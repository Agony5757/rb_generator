# RB Generator

[![PyPI version](https://badge.fury.io/py/rb-generator.svg)](https://badge.fury.io/py/rb-generator)
[![Test](https://github.com/Agony5757/rb_generator/actions/workflows/test.yml/badge.svg?branch=master)](https://github.com/Agony5757/rb_generator/actions/workflows/test.yml)

A high-performance C++ Python extension for quantum circuit Randomized Benchmarking (RB). It generates random Clifford sequences and their inverses for 2-qubit (RB22) and 4-qubit (RB44) systems.

## Features

- **High Performance**: C++17 core with OpenMP parallelism
- **2-Qubit RB (RB22)**: 11520 Clifford group elements, 7 basic generators
- **4-Qubit RB (RB44)**: Large Clifford group, 14 basic generators including CZ
- **Interleaved RB (IRB)**: Support for benchmarking specific gates
- **Auto Data Loading**: Automatic download or local generation of Clifford tables
- **Cross Platform**: Linux and Windows support

## Installation

### From PyPI (Recommended)

```bash
pip install rb_generator
```

### From Source

```bash
pip install .
```

### Editable Install

```bash
pip install -e .
```

## Requirements

- Python >= 3.10
- NumPy
- C++ compiler with C++17 support (g++ recommended on Linux)
- CMake >= 3.15
- OpenMP (optional, for parallel processing)

## Quick Start

```python
import rb_generator
import numpy as np

# Load Clifford tables (auto-downloads from GitHub if not present)
rb22, rb44 = rb_generator.load_table()

# Generate a random RB sequence
clifford_length = 100
cliffords = np.random.randint(0, rb22.N, clifford_length)

# Get operation sequence and its inverse
sequence, inverse_sequence = rb22.get_full_sequence_and_inverse_sequence(cliffords)

print(f"Sequence: {sequence}")
print(f"Inverse: {inverse_sequence}")

# Verify the sequence correctness
assert rb_generator.rb22_checker(sequence, inverse_sequence)
```

## API Reference

### Loading Tables

```python
rb_generator.load_table(path=".", from_="github")
```

Load rb22 and rb44 Clifford tables.

- **path**: Directory containing (or where to place) `rb22.dat` and `rb44.dat`
- **from_**: Data source when files are not present
  - `"github"`: Download precomputed tables from GitHub releases (default)
  - `"local"`: Generate tables locally using C++ generators

Returns a tuple `(RB22, RB44)` instances with loaded tables.

### RB22 (2-Qubit)

```python
rb22 = rb_generator.RB22()
rb22.load_from_file("./rb22.dat")

# Generate sequence
cliffords = np.random.randint(0, rb22.N, length)
sequence, inverse = rb22.get_full_sequence_and_inverse_sequence(cliffords)

# IRB: Get available special gates
special_gates = rb22.get_special_operations_str()
# Returns: ['I', 'X', 'Y', 'X/2', 'Y/2', '-X/2', '-Y/2']
```

### RB44 (4-Qubit)

```python
rb44 = rb_generator.RB44()
rb44.load_from_file("./rb44.dat")

# Generate sequence
cliffords = np.random.randint(0, rb44.N, length)
sequence, inverse = rb44.get_full_sequence_and_inverse_sequence(cliffords)

# IRB: Get available special gates
special_gates = rb44.get_special_operations_str()
# Returns: ['I', 'XI', 'IX', 'YI', 'IY', 'X/2_I', '-X/2_I', 'I_X/2',
#           'I_-X/2', 'Y/2_I', '-Y/2_I', 'I_Y/2', 'I_-Y/2', 'CZ']
```

### Interleaved RB (IRB)

```python
# Get special gate indices
special_names = rb22.get_special_operations_str()
special_indices = rb22.get_special_operations()

# Interleave a specific gate (e.g., 'X') between random Cliffords
interleaved_gate = special_indices[special_names.index('X')]
gate_sequences = np.reshape(
    np.vstack([cliffords, np.full(length, interleaved_gate)]),
    (-1,)
)
sequence, inverse = rb22.get_full_sequence_and_inverse_sequence(gate_sequences)
```

### Verification

```python
# Verify that sequence + inverse = identity
rb_generator.rb22_checker(sequence, inverse_sequence)
rb_generator.rb44_checker(sequence, inverse_sequence)
```

## Clifford Group Generators

### RB22 (2-qubit, 7 generators)

Defined in [clifford22.h](RBGeneratorCpp/src/clifford22.h):

| Index | Generator | Description |
|-------|-----------|-------------|
| 0 | I22 | Identity |
| 1 | X | Pauli-X |
| 2 | Y | Pauli-Y |
| 3 | SX | sqrt(X) |
| 4 | SY | sqrt(Y) |
| 5 | SXdag | sqrt(X) dagger |
| 6 | SYdag | sqrt(Y) dagger |

### RB44 (4-qubit, 14 generators)

Defined in [clifford44.h](RBGeneratorCpp/src/clifford44.h):

| Index | Generator | Description |
|-------|-----------|-------------|
| 0 | I44 | Identity |
| 1 | XI | X on qubit 0 |
| 2 | IX | X on qubit 1 |
| 3 | YI | Y on qubit 0 |
| 4 | IY | Y on qubit 1 |
| 5 | SX_I | sqrt(X) on qubit 0 |
| 6 | SXdag_I | sqrt(X) dagger on qubit 0 |
| 7 | I_SX | sqrt(X) on qubit 1 |
| 8 | I_SXdag | sqrt(X) dagger on qubit 1 |
| 9 | SY_I | sqrt(Y) on qubit 0 |
| 10 | SYdag_I | sqrt(Y) dagger on qubit 0 |
| 11 | I_SY | sqrt(Y) on qubit 1 |
| 12 | I_SYdag | sqrt(Y) dagger on qubit 1 |
| 13 | CZ | Controlled-Z |

## Data Files

The package requires precomputed Clifford group data files:
- `rb22.dat` (~500MB)
- `rb44.dat` (~500MB)

These files are not included in the PyPI wheel but are:
1. Automatically downloaded from GitHub releases on first use
2. Or can be generated locally (slower, requires C++ generator)

```python
# Force local generation
rb22, rb44 = rb_generator.load_table(from_="local")
```

## Building

```bash
# Build wheel
python -m build

# Run tests
pytest test/test.py -v
```

## License

Apache-2.0
