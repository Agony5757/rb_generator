import time
import rb_generator
import numpy as np

def test_rb22():
    rb22 = rb_generator.RB22()
    if not rb22.load_from_file('./rb22.dat'):
        print('load failed')
        exit(0)
    clifford_length = 100
    cliffords = np.random.randint(0, rb22.N, clifford_length)
    sequence, inverse_sequence = rb22.get_full_sequence_and_inverse_sequence(cliffords)
    if not rb_generator.rb22_checker(sequence, inverse_sequence):
        raise RuntimeError("Checker not passed.")
    print(sequence)
    print(inverse_sequence)

    
def test_rb44():

    rb44 = rb_generator.RB44()
    if not rb44.load_from_file('./rb44.dat'):
        print('load failed')
        exit(0)
    clifford_length = 100
    cliffords = np.random.randint(0, rb44.N, clifford_length)
    sequence, inverse_sequence = rb44.get_full_sequence_and_inverse_sequence(cliffords)
    if not rb_generator.rb44_checker(sequence, inverse_sequence):
            raise RuntimeError("Checker not passed.")
    print(sequence)
    print(inverse_sequence)
    

def benchmark_rb22(n_episodes = 1000, n_length = 1000):

    rb22 = rb_generator.RB22()
    if not rb22.load_from_file('./rb22.dat'):
        print('load failed')
        exit(0)

    t1 = time.time()
    for i in range(n_episodes):
        cliffords = np.random.randint(0, rb22.N, n_length)
        sequence, inverse_sequence = rb22.get_full_sequence_and_inverse_sequence(cliffords)
        if not rb_generator.rb22_checker(sequence, inverse_sequence):
            raise RuntimeError("Checker not passed.")
    t2 = time.time()

    print(f'Duration: {t2-t1} seconds, for {n_episodes} episodes (Clifford length = {n_length}).')

def benchmark_rb44(n_episodes = 1000, n_length = 1000):

    rb44 = rb_generator.RB44()
    if not rb44.load_from_file('./rb44.dat'):
        print('load failed')
        exit(0)

    t1 = time.time()
    for i in range(n_episodes):
        cliffords = np.random.randint(0, rb44.N, n_length)
        sequence, inverse_sequence = rb44.get_full_sequence_and_inverse_sequence(cliffords)
        if not rb_generator.rb44_checker(sequence, inverse_sequence):
            raise RuntimeError("Checker not passed.")
    t2 = time.time()

    print(f'Duration: {t2-t1} seconds, for {n_episodes} episodes (Clifford length = {n_length}).')

def test_irb22(interleaved_gate = 'X', n_episodes = 1000, n_length = 1000):

    rb22 = rb_generator.RB22()
    
    if not rb22.load_from_file('./rb22.dat'):
        print('load failed')
        exit(0)

    special_gate_names = rb22.get_special_operations_str()
    special_gates = rb22.get_special_operations()
    # print(special_gate_names)
    # print(special_gates)
    if interleaved_gate not in special_gate_names:
        raise RuntimeError('Invalid interleaved gate.')
    
    interleaved_gate_pos = special_gate_names.index(interleaved_gate)
    special_gate_names.index(interleaved_gate)
    interleaved_group_pos = special_gates[interleaved_gate_pos]

    t1 = time.time()
    for i in range(n_episodes):
        cliffords = np.random.randint(0, rb22.N, n_length)
        interleaved_gate = np.ones(n_length, dtype=np.int32) * interleaved_group_pos
        gate_sequences = np.reshape(np.vstack([cliffords, interleaved_gate]), (-1,))
        sequence, inverse_sequence = rb22.get_full_sequence_and_inverse_sequence(gate_sequences)
        if not rb_generator.rb22_checker(sequence, inverse_sequence):
            raise RuntimeError("Checker not passed.")
    t2 = time.time()

    print(f'Test passed. Duration: {t2-t1} seconds, for {n_episodes} episodes (Clifford length = {n_length}).')

def test_irb44(interleaved_gate = 'CZ', n_episodes = 1000, n_length = 1000):

    rb44 = rb_generator.RB44()
    
    if not rb44.load_from_file('./rb44.dat'):
        print('load failed')
        exit(0)

    special_gate_names = rb44.get_special_operations_str()
    special_gates = rb44.get_special_operations()
    # print(special_gate_names)
    # print(special_gates)
    if interleaved_gate not in special_gate_names:
        raise RuntimeError('Invalid interleaved gate.')
    
    interleaved_gate_pos = special_gate_names.index(interleaved_gate)
    special_gate_names.index(interleaved_gate)
    interleaved_group_pos = special_gates[interleaved_gate_pos]

    t1 = time.time()
    for i in range(n_episodes):
        cliffords = np.random.randint(0, rb44.N, n_length)
        interleaved_gate = np.ones(n_length, dtype=np.int32) * interleaved_group_pos
        gate_sequences = np.reshape(np.vstack([cliffords, interleaved_gate]).T, (-1,))
        sequence, inverse_sequence = rb44.get_full_sequence_and_inverse_sequence(gate_sequences)
        if not rb_generator.rb44_checker(sequence, inverse_sequence):
            raise RuntimeError("Checker not passed.")
    t2 = time.time()

    print(f'Test passed. Duration: {t2-t1} seconds, for {n_episodes} episodes (Clifford length = {n_length}).')

if __name__ == '__main__':
    # test_rb22()
    benchmark_rb22()
    test_irb22('I')
    test_irb22('X')
    test_irb22('Y')
    test_irb22('X/2')
    test_irb22('Y/2')
    test_irb22('-X/2')
    test_irb22('-Y/2')
    # test_rb44()
    benchmark_rb44()
    test_irb44('I')
    test_irb44('CZ')