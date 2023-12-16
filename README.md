# Randomized benchmarking for 1-qubit and 2-qubit 

## Install
```
pip install rb_generator
```

## Usage

See [test](test/test.py).

### rb22.get_full_sequence_and_inverse_sequence
```python
rb22 = rb_generator.RB22()

# Open file first
if not rb22.load_from_file('./rb22.dat'):
    print('load failed')
    exit(0)

# Define the clifford length
clifford_length = 100

# Generate a random sequence (each element in [0, N-1], N is the size of the clifford group)
cliffords = np.random.randint(0, rb22.N, clifford_length)
sequence, inverse_sequence = rb22.get_full_sequence_and_inverse_sequence(cliffords)

print(sequence)
print(inverse_sequence)
# You will obtain the following two outputs

# sequence = [3, 6, 5, 4, 1, 4, 4, 5, 2, 5, 3, 4, 3, 3, 6, 5, 6, 5, 4, 5, 1, 6, 3, 4, 5, 4, 4, 3, 5, 4, 3, 4, 5, 4, 6, 2, 3, 3, 4, 1, 4, 6, 5, 3, 3, 6, 3, 1, 6, 3, 6, 6, 5, 4, 3, 3, 4, 5, 4, 5, 5, 6, 4, 5, 0, 6, 5, 2, 5, 5, 4, 3, 6, 5, 2, 3, 5, 3, 6, 3, 4, 3, 4, 5, 3, 6, 6, 5, 3, 4, 5, 4, 3, 0, 3, 4, 3, 4, 2, 3, 2, 5, 1, 6, 4, 5, 6, 5, 3, 3, 6, 3, 5, 6, 5, 4, 5, 4, 3, 6, 3, 6, 5, 5, 4, 4, 5, 2, 2, 2, 1, 6, 6, 2, 3, 4, 3, 4, 5, 6, 6, 3, 6, 3, 4, 1, 4, 3, 6, 3, 0, 6, 5, 4, 5, 6, 3, 3, 4, 3, 4, 5, 6, 5, 6, 3, 2, 3, 6, 3, 0, 1, 6, 5, 6, 2, 3, 5, 4, 5, 6, 5, 6, 5]

# inverse_sequence = [5, 6]
```

"sequence" means the operation sequence for generating the RB sequence.
"inverse_sequence" means the operation sequence for inverting the sequence to identity.

### rb44.get_full_sequence_and_inverse_sequence
Same as the version of rb22.

## For the definition of basic generators
They are generated by C++ enums. Find them in [clifford22.h](RBGeneratorCpp/src/clifford22.h) and [clifford44.h](RBGeneratorCpp/src/clifford44.h)

```C++
enum Generator22Enum : int
{
    Generator_I22,
    Generator_X,
    Generator_Y,
    Generator_SX,
    Generator_SY,
    Generator_SXdag,
    Generator_SYdag,
};
```

```C++
enum Generator44Enum : int
{
    Generator_I44,
    Generator_XI,
    Generator_IX,
    Generator_YI,
    Generator_IY,
    Generator_SX_I,
    Generator_SXdag_I,
    Generator_I_SX,
    Generator_I_SXdag,
    Generator_SY_I,
    Generator_SYdag_I,
    Generator_I_SY,
    Generator_I_SYdag,
    Generator_CZ,
};
```
