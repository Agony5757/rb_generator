# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RB Generator is a C++ Python extension (pybind11) for quantum circuit randomized benchmarking. It generates random Clifford sequences and their inverses for 2-qubit (RB22) and 4-qubit (RB44) systems.

## Build Commands

```bash
# Install from source (builds C++ extension via CMake)
pip install .

# Editable install
pip install -e .

# Build wheel
python -m build

# Run tests (requires rb22.dat and rb44.dat data files in project root)
pytest test/test.py -v
```

## Versioning

Version is managed by `setuptools_scm` — derived automatically from git tags (`v0.0.5` → `0.0.5`). The generated `rb_generator/_version.py` is written at build time and exposes `__version__`. Never edit `_version.py` manually.

## Architecture

- **C++ core** (`RBGeneratorCpp/src/`): Clifford group implementations (`clifford22.h`, `clifford44.h`) and utilities (`utils.h`). Uses C++17, OpenMP for parallelism.
- **pybind11 bindings** (`RBGeneratorCpp/Pybinder/`): Exposes C++ classes to Python. `RBGeneratorCpp.cpp` is the main module entry point; `rb22.h` and `rb44.h` wrap the respective Clifford classes.
- **Vendored dependencies** (`RBGeneratorCpp/ThirdParty/`): pybind11 and fmt are bundled — do not add external dependencies on these.
- **Python package** (`rb_generator/__init__.py`): Imports the compiled extension and exposes `RB22`, `RB44`, `rb22_checker`, `rb44_checker`, `__version__`.
- **Build system**: `pyproject.toml` declares metadata and build dependencies. `setup.py` contains only the `CMakeExtension`/`CMakeBuild` classes for CMake integration. On Linux it forces `g++`/`gcc` compilers and prefers Ninja generator. Supports Python 3.10–3.14.

## Key Details

- Data files (`rb22.dat`, `rb44.dat`) contain precomputed Clifford group data (~500MB each), **not** in the repo or wheel. They are uploaded to GitHub releases for CI use. Runtime: first run auto-generates and caches locally.
- The Clifford group enums defining basic generators are in `clifford22.h` (2-qubit, 7 generators) and `clifford44.h` (4-qubit, 14 generators including CZ).
- `get_full_sequence_and_inverse_sequence` is the core method — takes a numpy array of Clifford indices and returns gate operation sequences with their inverses.
- Interleaved RB (IRB) support: `get_special_operations()` / `get_special_operations_str()` return gate mappings for interleaving specific gates between random Cliffords.
