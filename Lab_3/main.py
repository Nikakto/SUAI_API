from Lab_1.TransformationMatrix import TransformationMatrix
from Lab_3.dna import DNA

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

c_width = 5
y_height = 10


def _graph_pos(graph):

    pos = {}

    def _plot_neibs(graph, neibs, x, y):

        if not len(neibs):
            return

        line_width = sum(len(name) * c_width for name in neibs)
        name_width = line_width / len(neibs)
        for index, name in enumerate(neibs):

            name_x = x - line_width/2 + index * name_width
            pos[name] = (name_x, y * y_height)
            _plot_neibs(graph, list(graph.neighbors(name)), name_x, y-1)

    root = ''
    for node in graph.nodes:
        if 0 == len(list(graph.neighbors(node))):
            root = node
            break

    pos[root] = (0, 0)
    graph_rev = graph.reverse()
    _plot_neibs(graph_rev, list(graph_rev.neighbors(root)), 0, -1)

    return pos


def _print_matrix(matrix):

    for row in matrix:
        for value in row:
            s = '%.2f' % value
            print('%-8s' % s, end='')
        print()


if __name__ == '__main__':

    aligments = {}
    animals = []

    for animal_name, animal in DNA.items():
        print('\n\n', animal_name)
        animals.append(animal_name)

        for other_animal_name, other_animal in DNA.items():
            print('\tto', other_animal_name, ':')

            if animal_name == other_animal_name:
                continue
            elif other_animal_name in aligments.keys() and animal_name in aligments[other_animal_name]:
                continue

            total_aligment = 0
            for dna_code in other_animal.keys():
                print(dna_code)

                matrix = TransformationMatrix(animal[dna_code].replace('-', ''), other_animal[dna_code].replace('-', ''))
                matrix.balance_top_down(gap=1)
                aligment = matrix.get_aligment()

                total_aligment += aligment.count('S') + aligment.count('R') + aligment.count('I')

            if animal_name not in aligments.keys():
                aligments[animal_name] = {}
            aligments[animal_name][other_animal_name] = total_aligment

            if other_animal_name not in aligments.keys():
                aligments[other_animal_name] = {}
            aligments[other_animal_name][animal_name] = total_aligment

    matrix = np.zeros((len(animals), len(animals)), dtype=float)
    for animal_index, animal in enumerate(animals):
        for other_animal_index, other_animal in enumerate(animals):

            if animal_index == other_animal_index:
                matrix[animal_index][animal_index] = float('inf')
            else:
                matrix[animal_index][other_animal_index] = aligments[animal][other_animal]

    graph = nx.DiGraph()
    for animal in animals:
        graph.add_node(animal, animals=[animal, ])

    _print_matrix(matrix)
    while matrix.shape[0] > 1:

        index_0, index_1 = np.unravel_index(matrix.argmin(), matrix.shape)

        node_name = animals[index_0] + '__' + animals[index_1]
        node_animals = graph.nodes[animals[index_0]]['animals'] + graph.nodes[animals[index_1]]['animals']
        graph.add_node(node_name, animals=node_animals )

        edge_weight = matrix[index_0][index_1] / 2
        graph.add_edge(animals[index_0], node_name, weight=edge_weight)
        graph.add_edge(animals[index_1], node_name, weight=edge_weight)

        for column in range(matrix.shape[1]):

            if index_0 == column or index_1 == column:
                continue

            animal_x = animals[column]
            weight_node_0 = len(graph.nodes[animals[index_0]]['animals'])
            weight_node_1 = len(graph.nodes[animals[index_1]]['animals'])

            numerator = matrix[index_0][column] * weight_node_0 + matrix[index_1][column] * weight_node_1
            denumerator = weight_node_0 + weight_node_1
            new_dist = numerator / denumerator

            matrix[index_0][column] = new_dist
            matrix[column][index_0] = new_dist

        animals[index_0] = animals[index_0] + '__' + animals[index_1]
        animals.pop(index_1)

        matrix = np.delete(matrix, index_1, 0)
        matrix = np.delete(matrix, index_1, 1)

        print('\n\n')
        print(animals)
        _print_matrix(matrix)

    pos = _graph_pos(graph)
    nx.draw(graph, pos=pos, with_labels=True, node_color='w', font_size=10, linewidths=2)
    ax = plt.gca()  # to get the current axis

    x_min_name = min(pos, key=lambda x: pos[x][0])
    x_min = pos[x_min_name][0] - len(x_min_name) * c_width / 2

    x_max_name = max(pos, key=lambda x: pos[x][0])
    x_max = pos[x_max_name][0] + len(x_max_name) * c_width / 2

    ax.set_xlim(x_min, x_max)
    plt.show()