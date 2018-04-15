from Lab_2.stuff import _print_matrix, get_alignment, get_matrix

str_1 = 'TGTTACGG'
str_2 = 'GGTTGACTA'

use_file = True
file = './prince.txt'
str_1 = 'КОТОВАСЯ КИТОВ'
str_2 = 'ВЗРОСЛЫЙ'


min_match_length = 3

operations = {
    'S': -1000,
    'M': 1,
    'I': -1000,
    'R': -1000,
}

if __name__ == '__main__':

    aligment = ['', '']

    matrix = get_matrix(str_1, str_2, operations, debug=True)
    aligment[0], aligment[1], *_ = get_alignment(str_1, str_2, matrix)

    print('Best aligment:')
    for al in aligment:
        print(al)

    if use_file:
        with open(file, 'r') as f:
            str_1 = f.read().replace('\n', '').upper()
    matrix = get_matrix(str_1, str_2, operations, debug=False)

    match_length = min_match_length
    while matrix.max() > operations['M']:

        aligment[0], aligment[1], matrix_pos, match_length = get_alignment(str_1, str_2, matrix, remove_path=True)
        if match_length <= min_match_length:
            continue

        x_before = matrix_pos[1] - match_length - 10 if matrix_pos[1] - match_length - 10 > 0 else 0
        x_after = matrix_pos[1] + 10 if matrix_pos[1] + 10 < len(str_1) else len(str_1)

        underlined = "\033[4m%s\033[0m" % aligment[0]
        print('\nFound match in %d of str_1:' % (matrix_pos[1] - match_length + 1))
        print(
            str_1[x_before:(matrix_pos[1] - match_length + 1)] +
            underlined +
            str_1[matrix_pos[1]:x_after])

        print('Aligment:')
        for al in aligment:
            print(al)

    pass