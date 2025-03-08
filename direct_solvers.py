# Solve the problem!
from collections import deque
import copy

def bfs_shortest_distance(graph, start, target):
    queue = deque([start])
    distances = {start: 0}  # Distance from start to each node

    while queue:
        node = queue.popleft()

        if node[0] == target:  # Stop early if we reach the target
            return distances[node]

        for neighbor in graph.get(node, []):
            if neighbor not in distances:  # If not visited
                distances[neighbor] = distances[node] + 1
                queue.append(neighbor)

    return -1  # Return -1 if the target is unreachable


from itertools import product

def create_graph_from_strands(dna_list):
    assert dna_list

    from main import get_next_base, DNA_BASES


    # print(dna_list)

    # create vertices
    graph = {}
    len_of_strands = len(dna_list[0])
    num_of_strands = len(dna_list)
    # tuples = list(product(range(len_of_strands), repeat=len_of_strands))
    # tuples_as_lists = [list(tup) for tup in tuples]  # Convert tuples to list
    # tuples2 = list(product(tuples, DNA_BASES))
    # tuples_as_lists2 = [list(tup) for tup in tuples2]  # Convert tuples to list
    print("starting optimal solution")

    my_list = [list(tup) for tup in product(range(len_of_strands + 1), repeat=num_of_strands)]
    vertices_list = [list(tup) for tup in product(my_list, DNA_BASES)]
    print("wuff")

    #add edges
    for vertice in vertices_list:

        vertice_as_tuple = tuple(tuple(inner_list) for inner_list in vertice)
        graph[vertice_as_tuple] = []

        # print(vertice)
        list_of_possible_edges = []
        vertice_state = vertice[0]
        vertice_base = vertice[1]
        for i in range(num_of_strands):
            curr_strand = dna_list[i]
            curr_strand_next_base_index = vertice_state[i]
            if curr_strand_next_base_index >= len_of_strands:
                continue

            curr_strand_next_base = curr_strand[curr_strand_next_base_index]
            if curr_strand_next_base == vertice_base:
                # Add edge
                # Update the index
                vertice_to_connect_to_state = copy.deepcopy(vertice_state)
                vertice_to_connect_to_state[i] = vertice_to_connect_to_state[i] + 1
                # Update the base
                vertice_to_connect_to_base = get_next_base(vertice_base)

                vertice_to_connect_to = [vertice_to_connect_to_state, vertice_to_connect_to_base]
                vertice_to_connect_to_as_tupple = tuple(tuple(inner_list) for inner_list in vertice_to_connect_to)
                graph[vertice_as_tuple].append(vertice_to_connect_to_as_tupple)

        # Add an edge that does not move up in state, only to next base
        vertice_to_connect_to_same_state = copy.deepcopy(vertice_state)
        vertice_to_connect_to_same = [vertice_to_connect_to_same_state, get_next_base(vertice_base)]
        vertice_to_connect_to_same_state_as_tupple = tuple(tuple(inner_list) for inner_list in vertice_to_connect_to_same)
        graph[vertice_as_tuple].append(vertice_to_connect_to_same_state_as_tupple)

    # for key, value in graph.items():
    #     print(f"{key}: {value}")

    return graph

def get_optimal_row_solution(dna_list):
    assert dna_list
    from main import DNA_BASES

    graph = create_graph_from_strands(dna_list)
    starting_tupple = (tuple([0] * len(dna_list)), tuple(DNA_BASES[0]))
    target = tuple([len(dna_list[0])] * len(dna_list))
    return bfs_shortest_distance(graph,starting_tupple, target)