import matplotlib.pyplot as plt
import numpy as np

from Lab_1.TransformationMatrix import TransformationMatrix
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

    transformation = TransformationMatrix(DNA_0, DNA_1)

    # transformation.balance_bottom_up()
    # print('\n', transformation)
    # print('\nOperations [Recursive]: ', transformation.operations, end='\n\n')

    # transformation.balance_bottom_up(gap=-1, similarity_matrix=SIMILARITY_MATRIX)
    # transformation.balance_bottom_up(gap=-1, similarity_matrix=transformation.generate_similarity_matrix())
    # print('PATH: ', transformation.get_path(maximum=True))
    # transformation.balance_bottom_up()
    # print('PATH: ', transformation.get_path())
    # print('\n', transformation)
    # print('\nOp/erations [bottom-up]: ', transformation.operations, end='\n\n')

    # transformation.balance_top_down(gap=-1, similarity_matrix=SIMILARITY_MATRIX)
    # transformation.balance_top_down(gap=-1, similarity_matrix=transformation.generate_similarity_matrix())
    # print('PATH: ', transformation.get_path(maximum=True))
    # print('ALIGMENT:', transformation.get_aligment(maximum=True))
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