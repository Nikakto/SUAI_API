from collections import UserList

from Lab_1.const import NUCLEOTIDES


class TransformationMatrix(UserList):

    DNA = None
    operations = 0

    def __init__(self, DNA_0, DNA_1, similarity_matrix=None):
        super().__init__()

        DNA_0 = ' ' + DNA_0
        DNA_1 = ' ' + DNA_1
        self.DNA = [DNA_0, DNA_1]

    def __str__(self):

        empty_space = [' ' * 5]

        header = ['%-5s' % f'{char}' for char in self.DNA[0]]
        string_rows = [''.join(empty_space + header)]

        for row_index, row in enumerate(self):

            string_row = ['%-5s' % f'{self.DNA[1][row_index]}']
            for column in row:

                if column >= 0:
                    value = ' %-4s' % f'{column}'
                else:
                    value = '%-5s' % f'{column}'
                string_row.append(value)

            string_rows.append(''.join(string_row))

        return '\n'.join(string_rows)

    def _diff(self, row, column, similarity_matrix=None):

        if similarity_matrix:
            nucleotide_column = self.DNA[0][column]
            nucleotide_row = self.DNA[1][row]
            return similarity_matrix[nucleotide_row][nucleotide_column]
        else:
            return self.DNA[0][column] != self.DNA[1][row]

    def _get_minimum(self, row, column, maximum):

        can_insert = column - 1 >= 0
        can_remove = row - 1 >= 0

        if not can_insert and not can_remove:
            return None

        cells = []
        if can_insert:
            cells.append((row, column - 1))

        if can_remove:
            cells.append((row - 1, column))

        if can_insert and can_remove:
            cells.append((row - 1, column - 1))

        if maximum:
            return max(cells, key=lambda cell: (self[cell[0]][cell[1]], cell[0], cell[1]))
        else:
            return min(cells, key=lambda cell: (self[cell[0]][cell[1]], cell[0], cell[1]))

    def _get_minimum_price(self, row, column, gap=1, similarity_matrix=None):

        can_insert = column - 1 >= 0
        can_remove = row - 1 >= 0

        if not can_insert and not can_remove:
            return 0

        values = []

        if can_insert:
            values.append(self[row][column - 1] + gap)

        if can_remove:
            values.append(self[row - 1][column] + gap)

        if can_insert and can_remove:
            diff = self._diff(row, column, similarity_matrix)
            values.append(self[row - 1][column - 1] + diff)

        if similarity_matrix:
            # look max similarity
            return max(values)
        else:
            # look minimum operation cost
            return min(values)

    def balance_bottom_up(self, row=None, column=None, gap=1, similarity_matrix=None):

        self.operations += 1
        if column is None and row is None:

            self.clear()
            self.operations = 0

            for i in range(len(self.DNA[1])):
                self.append([float('inf')] * len(self.DNA[0]))

            if similarity_matrix is not None:

                for column in range(len(self.DNA[0])):
                    self[0][column] = gap * column

                for row in range(len(self.DNA[1])):
                    self[row][0] = gap * row

            row = len(self.DNA[1]) - 1
            column = len(self.DNA[0]) - 1

        if self[row][column] == float('inf'):

            if column == 0:
                self[row][column] = row
            elif row == 0:
                self[row][column] = column
            else:

                insert = self.balance_bottom_up(row, column - 1, gap, similarity_matrix) + gap
                delete = self.balance_bottom_up(row - 1, column, gap, similarity_matrix) + gap

                diff = self._diff(row, column, similarity_matrix)
                substitute = self.balance_bottom_up(row - 1, column - 1, gap, similarity_matrix) + diff

                if similarity_matrix:
                    # look max similarity
                    self[row][column] = max(insert, delete, substitute)
                else:
                    # look minimum operation cost
                    self[row][column] = min(insert, delete, substitute)

        return self[row][column]

    def balance_top_down(self, gap=1, similarity_matrix=None):

        self.clear()
        self.operations = 0

        for i in range(len(self.DNA[1])):
            self.append([float('inf')] * len(self.DNA[0]))

        row = 0
        column = 0
        self[row][column] = 0

        def get_diagonal(row_d, column_d):
            result = []
            while row_d < len(self.DNA[1]) and column_d >= 0:
                result.append((row_d, column_d))
                row_d += 1
                column_d -= 1

            return result

        row, column = 0, 0
        while column < len(self.DNA[0]) and row < len(self.DNA[1]):

            cells = get_diagonal(row, column)
            for cell_row, cell_column in cells:
                self.operations += 1
                self[cell_row][cell_column] = self._get_minimum_price(cell_row, cell_column, gap, similarity_matrix)

            if column < len(self.DNA[0]) - 1:
                column += 1
            else:
                row += 1

    def generate_similarity_matrix(self):

        similarity_matrix = {

            nucleotide_column: {
                nucleotide_row: 1 - 2*int(nucleotide_row != nucleotide_column) for nucleotide_row in set(self.DNA[0])
            }
            for nucleotide_column in set(self.DNA[1])
        }

        return similarity_matrix

    def get_aligment(self, maximum=False):

        path = self.get_path(maximum=maximum)
        aligment = ''
        for i in range(len(path)-1):
            cell_from, cell_to = path[i], path[i+1]

            if cell_from[0] == cell_to[0]:
                aligment += 'I'
            elif cell_from[1] == cell_to[1]:
                aligment += 'D'
            elif self.DNA[1][cell_from[0]] == self.DNA[0][cell_from[1]]:
                aligment += 'M'
            else:
                aligment += 'R'

        return ''.join(reversed(aligment))

    def get_path(self, row=None, column=None, maximum=False):

        row = len(self.DNA[1]) - 1 if row is None else row
        column = len(self.DNA[0]) - 1 if column is None else column

        previous = self._get_minimum(row, column, maximum)
        if previous is None:
            return [(row, column)]
        else:
            return [(row, column)] + self.get_path(previous[0], previous[1], maximum)

    def print_aligment(self, maximum=False):

        print(f'FROM: {self.DNA[1]}')
        print(f'TO  : {self.DNA[0]}')
        print(f'{"="*20}')

        aligment = self.get_aligment(maximum=maximum)