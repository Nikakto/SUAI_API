import matplotlib.pyplot as plt
import numpy as np

from Lab_1.TransformationMatrix import TransformationMatrix, hirschberg
from Lab_1.const import NUCLEOTIDES, SIMILARITY_MATRIX


if __name__ == '__main__':

    DNA_0 = 'GAGGTAGGGAAC'
    DNA_1 = 'TGACCCCACG'

    # DNA_0 = ''.join(np.random.choice(NUCLEOTIDES, 12))
    # DNA_1 = ''.join(np.random.choice(NUCLEOTIDES, 10))

    # debug
    DNA_0 = 'EDITING'
    DNA_1 = 'DISTANCE'

    # debug2
    # DNA_0 = 'GTTAC'
    # DNA_1 = 'GACGT'

    # debug3
    # DNA_0 = 'GCATGCU'
    # DNA_1 = 'GATTACA'

    # debug4
    DNA_0 = 'AGTACGCA'
    DNA_1 = 'TATGC'

    transformation = TransformationMatrix(DNA_0, DNA_1)

    # transformation.balance_bottom_up()
    # print('\n', transformation)
    # print('\nOperations [Recursive]: ', transformation.operations, end='\n\n')

    # transformation.balance_bottom_up(gap=-1, similarity_matrix=SIMILARITY_MATRIX)
    # transformation.balance_bottom_up(gap=-1, similarity_matrix=transformation.generate_similarity_matrix())
    # print('PATH: ', transformation.get_path(maximum=True))
    transformation.balance_bottom_up()
    print('PATH: ', transformation.get_path())
    print('\n', transformation)
    print('\nOp/erations [bottom-up]: ', transformation.operations, end='\n\n')

    # transformation.balance_top_down(gap=-1, similarity_matrix=SIMILARITY_MATRIX)
    transformation.balance_top_down(gap=-1, similarity_matrix=transformation.generate_similarity_matrix())
    print('PATH: ', transformation.get_path(maximum=True))
    print('ALIGMENT:', transformation.get_aligment(maximum=True))


    # THIS UNCOM

    transformation.balance_top_down()
    print('PATH: ', transformation.get_path())
    print('ALIGMENT:', transformation.get_aligment())
    print('\n', transformation)
    print('\nOperations [top-down]: ', transformation.operations, end='\n\n')
    print('Done')

    recurs = [0]
    linear = [0]
    for n in range(1, 25):

        DNA_0 = ''.join(np.random.choice(NUCLEOTIDES, n))
        DNA_1 = ''.join(np.random.choice(NUCLEOTIDES, n))

        transformation = TransformationMatrix(DNA_0, DNA_1)

        transformation.balance_bottom_up()
        recurs.append(transformation.operations)

        transformation.balance_top_down()
        linear.append(transformation.operations)

    line_top_down, = plt.plot(range(25), linear, 'k-', label='top_down')
    line_bottom_up, = plt.plot(range(25), recurs, 'k--', label='bottom_up')
    plt.legend(handles=[line_top_down, line_bottom_up])
    plt.xlabel('n, n=m')
    plt.ylabel('operations')
    plt.show()


    # DOP
    l0, l1 = hirschberg(DNA_0, DNA_1)
    print('%s\n%s' % (l0, l1))

    x = []
    y0 = []
    yh = []

    t0 = []
    th = []

    for n in range(5, 25):

        import random
        from time import perf_counter
        DNA_0 = ''.join(np.random.choice(NUCLEOTIDES, n))
        DNA_1 = list(DNA_0)
        for i in range(len(DNA_0)):
            if random.random() < 0.5:
                DNA_1[i] = random.choice(NUCLEOTIDES)
        DNA_1 = ''.join(DNA_1)

        print('\n\nTEST %d' % n)
        t = perf_counter()
        print('%s\n%s' % hirschberg(DNA_0, DNA_1, sub_=2))
        th.append(perf_counter() - t)

        matrix = TransformationMatrix(DNA_0, DNA_1)
        t = perf_counter()
        matrix.balance_top_down(gap=1)
        t0.append(perf_counter() - t)
        matrix.get_aligment()

        import sys
        x.append(n)
        y0.append(72+32+8*sum(len(row) for row in matrix))
        yh.append(n*2*24)

    top_down, = plt.plot(x, t0, 'k-', label='top_down')
    hirschberg, = plt.plot(x, th, 'k--', label='Hirschberg')
    plt.legend(handles=[top_down, hirschberg])
    plt.xlabel('n, n=m')
    plt.ylabel('time (sec)')
    plt.show()

    top_down, = plt.plot(x, y0, 'k-', label='top_down')
    hirschberg, = plt.plot(x, yh, 'k--', label='Hirschberg')
    plt.legend(handles=[top_down, hirschberg])
    plt.xlabel('n, n=m')
    plt.ylabel('memory (byte)')
    plt.show()