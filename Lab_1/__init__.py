import numpy as np


from Lab_1.TransformationMatrix import TransformationMatrix
from Lab_1.const import NUCLEOTIDES, SIMILARITY_MATRIX


if __name__ == '__main__':

    n, m = 10, 10
    DNA_0 = np.random.choice(NUCLEOTIDES, n)
    DNA_1 = np.random.choice(NUCLEOTIDES, m)

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

    transformation.balance_bottom_up(gap=-1, similarity_matrix=transformation.generate_similarity_matrix())
    # transformation.balance_recursive(gap=-1, similarity_matrix=SIMILARITY_MATRIX)
    print('\n', transformation)
    print('\nOperations [bottom-up]: ', transformation.operations, end='\n\n')

    transformation.balance_top_down(gap=-1, similarity_matrix=transformation.generate_similarity_matrix())
    # transformation.balance_recursive(gap=-1, similarity_matrix=SIMILARITY_MATRIX)
    print('\n', transformation)
    print('\nOperations [top-down]: ', transformation.operations, end='\n\n')
    print('Done')
