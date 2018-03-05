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

    def _diff(self, row, column, similarity_matrix=None):

        if similarity_matrix:
            nucleotide_column = self.DNA[0][column]
            nucleotide_row = self.DNA[1][row]
            return similarity_matrix[nucleotide_row][nucleotide_column]
        else:
            return self.DNA[0][column] != self.DNA[1][row]

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
                delete = self.balance_bottom_up(row -1, column, gap, similarity_matrix) + gap

                diff = self._diff(row, column, similarity_matrix)
                substitute = self.balance_bottom_up(row - 1, column - 1, gap, similarity_matrix) + diff

                if similarity_matrix:
                    # look max similarity
                    self[row][column] = max(insert, delete, substitute)
                else:
                    # look minimum operation cost
                    self[row][column] = min(insert, delete, substitute)

        return self[row][column]

    def balance_top_down(self, row=None, column=None, gap=1, similarity_matrix=None):

        self.operations += 1

        if column is None and row is None:

            self.clear()
            self.operations = 0

            for i in range(len(self.DNA[1])):
                self.append([float('inf')] * len(self.DNA[0]))

            # if similarity_matrix is not None:
            #
            #     for column in range(len(self.DNA[0])):
            #         self[0][column] = gap * column
            #
            #     for row in range(len(self.DNA[1])):
            #         self[row][0] = gap * row

            row = 0
            column = 0
            self[row][column] = 0

        # print(self)
        can_insert = column+1 < len(self[row])
        can_remove = row+1 < len(self)
        insert = self[row][column] + gap
        delete = self[row][column] + gap

        if can_insert and can_remove:
            diff = self._diff(row + 1, column + 1, similarity_matrix)
            substitute = self[row][column] + diff
            if self[row + 1][column + 1] == float('inf') or self[row + 1][column + 1] < substitute:
                self[row + 1][column + 1] = substitute
                self.balance_top_down(row + 1, column + 1, gap, similarity_matrix)

        if can_remove and (self[row + 1][column] == float('inf') or self[row + 1][column] < delete):
            self[row + 1][column] = delete
            self.balance_top_down(row + 1, column, gap, similarity_matrix)

        if can_insert and (self[row][column + 1] == float('inf') or self[row][column + 1] < insert):
            self[row][column + 1] = insert
            self.balance_top_down(row, column + 1, gap, similarity_matrix)

    def generate_similarity_matrix(self):

        similarity_matrix = {

            nucleotide_column: {
                nucleotide_row: 1 - 2*int(nucleotide_row != nucleotide_column) for nucleotide_row in set(self.DNA[0])
            }
            for nucleotide_column in set(self.DNA[1])
        }

        return similarity_matrix

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