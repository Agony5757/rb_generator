import time
import rb_generator
import numpy as np

def testrb22():
    rb22 = rb_generator.RB22()
    if not rb22.load_from_file('./rb22.dat'):
        print('load failed')
        exit(0)
    clifford_length = 100
    cliffords = np.random.randint(0, rb22.N, clifford_length)
    sequence, inverse_sequence = rb22.get_full_sequence_and_inverse_sequence(cliffords)

    print(sequence)
    print(inverse_sequence)

    
def testrb44():

    rb44 = rb_generator.RB44()
    if not rb44.load_from_file('./rb44.dat'):
        print('load failed')
        exit(0)
    clifford_length = 100
    cliffords = np.random.randint(0, rb44.N, clifford_length)
    sequence, inverse_sequence = rb44.get_full_sequence_and_inverse_sequence(cliffords)

    print(sequence)
    print(inverse_sequence)

def benchmark_rb44(n_episodes = 1000, n_length = 1000):

    rb44 = rb_generator.RB44()
    if not rb44.load_from_file('./rb44.dat'):
        print('load failed')
        exit(0)

    t1 = time.time()
    for i in range(n_episodes):
        cliffords = np.random.randint(0, rb44.N, n_length)
        sequence, inverse_sequence = rb44.get_full_sequence_and_inverse_sequence(cliffords)

    t2 = time.time()

    print(f'Duration: {t2-t1} seconds, for {n_episodes} episodes (Clifford length = {n_length}).')

if __name__ == '__main__':
    benchmark_rb44()