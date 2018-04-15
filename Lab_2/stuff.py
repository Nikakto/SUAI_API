import numpy as np
import scipy
import scipy.sparse as sp


def _get_back_path(matrix, row, column):

    if matrix[row-1, column-1] <= 0:
        return [(row, column), (row-1, column-1)]

    matrix_part = matrix[row-1:row+1, column-1:column+1]
    matrix_part[1, 1] = 0
    same_row, same_column = scipy.unravel_index(matrix_part.argmax(), matrix_part.shape)

    if same_row and same_column:
        return [(row, column)]

    point_row = row if same_row else row - 1
    point_column = column if same_column else column - 1

    return [(row, column), ] + _get_back_path(matrix, point_row, point_column)


def _print_matrix(str_1, str_2, matrix):

    # HEADER
    print(' ' * 4 * 2, end='')
    print(''.join(["%-4s" % c for c in str_1]))

    width = matrix.shape[1]

    # FIRST LINE
    print(' ' * 4, end='')
    print(''.join(
        ["%-4d" % matrix[0, i] for i in range(width)]
    ))

    # OTHER LINE
    for index, c in enumerate(str_2):
        print('%-4s' % c, end='')
        print(''.join(
            ["%-4d" % matrix[index+1, i] for i in range(width)]
        ))


def get_alignment(str_1, str_2, matrix, remove_path=False):

    max_row, max_column = scipy.unravel_index(matrix.argmax(), matrix.shape)
    path = _get_back_path(matrix, max_row, max_column)

    aligment_1 = ''
    aligment_2 = ''

    for point_index in range(1, len(path)):

        row_0, column_0 = path[point_index]
        row_1, column_1 = path[point_index - 1]

        if row_0 != row_1 and column_0 == column_1:
            aligment_1 += '-'
            aligment_2 += str_2[row_0]
        elif row_0 == row_1 and column_0 != column_1:
            aligment_1 += str_1[column_0]
            aligment_2 += '-'
        else:
            aligment_1 += str_1[column_0]
            aligment_2 += str_2[row_0]

    if remove_path:

        x0, y0 = path[0]
        x1, y1 = path[-2]

        for i in range(x1, x0 + 1):
            matrix[i, y1:y0+1] = 0

        matrix.eliminate_zeros()

    return aligment_1[::-1], aligment_2[::-1], path[0], len(path)


def get_matrix(str_1, str_2, operations, debug=False):
    '''
    :param str_1: string destination
    :param str_2: string align to
    :return: str_2 to str_1
    '''

    shape = ( len(str_2) + 1, len(str_1) + 1, )
    matrix = sp.csr_matrix(shape)

    for row_index in range(1, matrix.shape[0]):
        for column_index in range(1, matrix.shape[1]):

            cost_sub = operations['S'] if str_1[column_index - 1] != str_2[row_index - 1] else operations['M']
            sub = matrix[row_index - 1, column_index - 1] + cost_sub

            ins = matrix[row_index, column_index - 1] + operations['I']
            rem = matrix[row_index - 1, column_index] + operations['R']

            matrix[row_index, column_index] = max(sub, ins, rem, 0)

    if debug:
        _print_matrix(str_1, str_2, matrix)

    return matrix