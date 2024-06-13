import xml.etree.ElementTree as ET
import networkx as nx
import random
import sys

class DNASequence:
    def __init__(self, length, elements_WS, elements_RY, start):
        self.length = length
        self.elements_WS = elements_WS
        self.elements_RY = elements_RY
        self.start = start

def load_from_xml(file_path):
    try:
        print("Attempting to load XML file from:", file_path)
        tree = ET.parse(file_path)
    except Exception as e:
        return f"Error loading XML: {str(e)} ({e.__class__.__name__})"
    root = tree.getroot()

    length = int(root.get('length'))
    start = root.get('start')

    elements_WS = []
    elements_RY = []
    for probe in root.findall('.//probe'):
        pattern = probe.get('pattern')
        for cell in probe.findall('cell'):
            if cell.text:
                if pattern.startswith('Z'):
                    elements_WS.append(cell.text)
                elif pattern.startswith('P'):
                    elements_RY.append(cell.text)
    return DNASequence(length, elements_WS, elements_RY, start)

def find_common_suffix(word1, word2):
    common_suffix = ""
    for i in range(1, min(len(word1), len(word2)) + 1):
        if word1[-i:] == word2[:i]:
            common_suffix = word1[-i:]
            break
    return common_suffix

def degenerate_sequence(sequence, degenerate_type):
    degenerate_map_WS = {'A': 'W', 'T': 'W', 'C': 'S', 'G': 'S'}
    degenerate_map_RY = {'A': 'R', 'G': 'R', 'C': 'Y', 'T': 'Y'}
    
    if degenerate_type == "WS":
        return ''.join(degenerate_map_WS.get(nuc, nuc) for nuc in sequence)
    elif degenerate_type == "RY":
        return ''.join(degenerate_map_RY.get(nuc, nuc) for nuc in sequence)
    else:
        return sequence

def merge_degenerate_nucleotides(nucleotide1, nucleotide2):
    # Inicjujemy pusty wynikowy nukleotyd
    merged_nucleotide = ""

    # Iterujemy po każdej literze z pierwszego nukleotydu
    for i in range(len(nucleotide1)):
        # Jeśli litera jest niezdegenerowana, dodaj ją do wyniku
        if nucleotide1[i] not in ('W', 'S'):
            merged_nucleotide += nucleotide1[i]
        # Jeśli litera jest zdegenerowana, sprawdź komplementarną część z drugiego nukleotydu
        elif nucleotide1[i] == 'W':
            if nucleotide2[i] == 'R':
                merged_nucleotide += 'A'
            elif nucleotide2[i] == 'Y':
                merged_nucleotide += 'T'
        elif nucleotide1[i] == 'S':
            if nucleotide2[i] == 'R':
                merged_nucleotide += 'G'
            elif nucleotide2[i] == 'Y':
                merged_nucleotide += 'C'

    return merged_nucleotide


def create_graph(sequences, start_sequence, degenerate_type):
    G = nx.DiGraph()

    G.add_node(start_sequence, label=start_sequence)

    for word1 in sequences:
        if word1 != start_sequence:
            G.add_node(word1, label=word1)
            common_suffix = find_common_suffix(degenerate_sequence(start_sequence, degenerate_type), degenerate_sequence(word1, degenerate_type))
            weight = len(common_suffix)
            if weight > 0:
                G.add_edge(start_sequence, word1, weight=weight)

    for word1 in sequences:
        for word2 in sequences:
            if word1 == word2 or word1 == start_sequence or word2 == start_sequence:
                continue
            res = find_common_suffix(degenerate_sequence(word1, degenerate_type), degenerate_sequence(word2, degenerate_type))
            weight = len(res)
            if weight > 0:
                G.add_edge(word1, word2, weight=weight)

    return G

def validate_degenerate_sequence(sequence, new_sequence, common_suffix_length, other_graph, other_deg_type):
    required_suffix = new_sequence[-1]
    required_prefix = degenerate_sequence(sequence[-common_suffix_length:], other_deg_type)

    for node in other_graph.nodes():
        node_label = other_graph.nodes[node]['label']
        if node_label[-1] == required_suffix and node_label.startswith(required_prefix) and not node_label.startswith(degenerate_sequence('A' * common_suffix_length, other_deg_type)):
            return True
    return False

def best_out_degree_neighbor_search(graph, current_node, max_length, current_sequence, errors, visited, other_graph, other_deg_type, best_solution):
    # Debug print to track the current state
    # print(f"Best out-degree neighbor search at node {current_node}, sequence: {current_sequence}, errors: {errors}")

    if len(current_sequence) >= max_length:
        if errors < best_solution[1]:
            best_solution[0] = current_sequence
            best_solution[1] = errors
        return

    # Znajdź sąsiadujące sekwencje
    neighbors = list(graph.neighbors(current_node))

    # Wybierz sąsiada na podstawie ilości wychodzących łuków
    neighbors_with_degrees = [(neighbor, graph.out_degree(neighbor)) for neighbor in neighbors]
    neighbors_with_degrees.sort(key=lambda x: x[1], reverse=True)

    # Przetasuj sąsiadów, aby uniknąć zawsze wybierania tych samych
    random.shuffle(neighbors_with_degrees)

    # Sprawdź każdego sąsiada w losowej kolejności
    for neighbor, _ in neighbors_with_degrees:
        edge_data = graph.get_edge_data(current_node, neighbor)
        weight = edge_data['weight']
        new_sequence = current_sequence + graph.nodes[neighbor]['label'][weight:]
        new_errors = errors + (len(graph.nodes[neighbor]['label']) - weight - 1)

        # Sprawdź, czy wierzchołek jest odwiedzony
        if neighbor in visited:
            continue

        # Sprawdź, czy istnieje dopełniający nukleotyd w drugim grafie
        if validate_degenerate_sequence(current_sequence, graph.nodes[neighbor]['label'], weight, other_graph, other_deg_type):
            # Znajdź dopasowany nukleotyd w drugim grafie
            complement_node = None
            for node in other_graph.nodes():
                node_label = other_graph.nodes[node]['label']
                if node_label.endswith(graph.nodes[neighbor]['label'][-1]):
                    complement_node = node
                    break
            if complement_node is not None:
                # Dodaj tylko końcówkę, która nie jest pokryta
                complement_label = other_graph.nodes[complement_node]['label']
                suffix_length = len(current_sequence)
                new_suffix = complement_label[suffix_length:]
                # print(f"Considered cell: {graph.nodes[neighbor]['label']}")
                # print(f"Found complementary cell from the other half of the graph: {complement_label}")
                # print(f"Added nucleotide: {new_suffix}")

                # Połącz obie komórki w jeden niezdegenerowany nukleotyd
                merged_nucleotide = merge_degenerate_nucleotides(graph.nodes[neighbor]['label'], complement_label)
                # print(f"Merged nucleotide: {merged_nucleotide}")
                merged_res = merged_nucleotide[weight:]
                # print(f"Should get added: {merged_nucleotide}")
                new_sequence = current_sequence + merged_res

                visited.add(neighbor)  # Oznacz wierzchołek jako odwiedzony
                best_out_degree_neighbor_search(graph, neighbor, max_length, new_sequence, new_errors, visited, other_graph, other_deg_type, best_solution)
                visited.remove(neighbor)  # Usuń oznaczenie wierzchołka jako odwiedzony

                # Dodaj niezdegenerowaną wersję nukleotydu dopiero po sprawdzeniu warunku validate_degenerate_sequence
                if len(best_solution[0]) >= max_length:
                    return



def dfs(graph, current_node, max_length, current_sequence, errors, visited, other_graph, other_deg_type, best_solution):
    # Debug print to track the current state
    # print(f"DFS at node {current_node}, sequence: {current_sequence}, errors: {errors}")

    if len(current_sequence) >= max_length:
        if errors < best_solution[1]:
            best_solution[0] = current_sequence
            best_solution[1] = errors
        return

    # Znajdź sąsiadujące sekwencje
    neighbors = list(graph.neighbors(current_node))

    # Sprawdź każdego sąsiada
    for neighbor in neighbors:
        edge_data = graph.get_edge_data(current_node, neighbor)
        weight = edge_data['weight']
        new_sequence = current_sequence + graph.nodes[neighbor]['label'][weight:]
        new_errors = errors + (len(graph.nodes[neighbor]['label']) - weight - 1)

        # Sprawdź, czy wierzchołek jest odwiedzony
        if neighbor in visited:
            continue

        visited.add(neighbor)  # Oznacz wierzchołek jako odwiedzony

        # Sprawdź, czy istnieje dopełniający nukleotyd w drugim grafie
        if validate_degenerate_sequence(current_sequence, graph.nodes[neighbor]['label'], weight, other_graph, other_deg_type):
            # Znajdź dopasowany nukleotyd w drugim grafie
            complement_node = None
            for node in other_graph.nodes():
                node_label = other_graph.nodes[node]['label']
                if node_label.endswith(graph.nodes[neighbor]['label'][-1]):
                    complement_node = node
                    break
            if complement_node is not None:
                # Dodaj tylko końcówkę, która nie jest pokryta
                complement_label = other_graph.nodes[complement_node]['label']
                suffix_length = len(current_sequence)
                new_suffix = complement_label[suffix_length:]

                # Połącz obie komórki w jeden niezdegenerowany nukleotyd
                merged_nucleotide = merge_degenerate_nucleotides(graph.nodes[neighbor]['label'], complement_label)
                merged_res = merged_nucleotide[weight:]
                new_sequence = current_sequence + merged_res

                dfs(graph, neighbor, max_length, new_sequence, new_errors, visited, other_graph, other_deg_type, best_solution)

        visited.remove(neighbor)  # Usuń oznaczenie wierzchołka jako odwiedzony

def solve(xml_file):
    try:
        probe = load_from_xml(xml_file)
        if isinstance(probe, str):
            return probe

        if not probe.elements_WS or not isinstance(probe.elements_WS[0], str):
            return "Error solving: Invalid DNA sequences in elements_WS"
        if not probe.elements_RY or not isinstance(probe.elements_RY[0], str):
            return "Error solving: Invalid DNA sequences in elements_RY"

        start_sequence = probe.start

        graph_WS = create_graph(probe.elements_WS, start_sequence, "WS")
        graph_RY = create_graph(probe.elements_RY, start_sequence, "RY")
        print(f"Start: {probe.start}, Length: {probe.length}")

        best_solution_WS = ["", float('inf')]
        best_solution_RY = ["", float('inf')]

        # Iteruj po wszystkich wierzchołkach w graph_WS
        for start_node_WS in graph_WS.nodes():
            dfs(graph_WS, start_node_WS, probe.length, start_sequence, 0, set([start_sequence]), graph_RY, "RY", best_solution_WS)

        # Iteruj po wszystkich wierzchołkach w graph_RY
        for start_node_RY in graph_RY.nodes():
            dfs(graph_RY, start_node_RY, probe.length, start_sequence, 0, set([start_sequence]), graph_WS, "WS", best_solution_RY)

        # Porównaj najlepsze rozwiązania z obu grafów i wybierz lepsze
        best_solution = min(best_solution_WS, best_solution_RY, key=lambda x: x[1])

        return best_solution, best_solution[1]  # Zwracamy najlepsze rozwiązanie i liczbę błędów
    except Exception as e:
        return f"Error solving: {str(e)} ({e.__class__.__name__})"


def heuristic_solve(xml_file):
    try:
        probe = load_from_xml(xml_file)
        if isinstance(probe, str):
            return probe

        if not probe.elements_WS or not isinstance(probe.elements_WS[0], str):
            return "Error solving: Invalid DNA sequences in elements_WS"
        if not probe.elements_RY or not isinstance(probe.elements_RY[0], str):
            return "Error solving: Invalid DNA sequences in elements_RY"

        start_sequence = probe.start

        graph_WS = create_graph(probe.elements_WS, start_sequence, "WS")
        graph_RY = create_graph(probe.elements_RY, start_sequence, "RY")
        print(f"Start: {probe.start}, Length: {probe.length}")

        best_solution_WS = ["", float('inf')]
        start_nodes_WS = [v for v in graph_WS.nodes(data=True) if v[1]['label'] == start_sequence]
        if not start_nodes_WS:
            raise ValueError(f"No node found with label {start_sequence} in graph_WS")
        for start_node_WS in start_nodes_WS:
            best_out_degree_neighbor_search(graph_WS, start_node_WS[0], probe.length, start_sequence, 0, set([start_sequence]), graph_RY, "RY", best_solution_WS)

        return best_solution_WS, best_solution_WS[1]  # Zwracamy najlepsze rozwiązanie i liczbę błędów
    except Exception as e:
        return f"Error solving: {str(e)} ({e.__class__.__name__})"


# Przykład wywołania
xml_file = sys.argv[1] if len(sys.argv) > 1 else 'input.xml'
mode = sys.argv[2] if len(sys.argv) > 2 else 'normal'

if mode == 'heuristic':
    result = heuristic_solve(xml_file)
else:
    result = solve(xml_file)

if isinstance(result, str):
    print(result)  # Wypisz komunikat o błędzie, jeśli wystąpił
else:
    best_solution, errors = result
    print("Best solution:", best_solution[0])
    print("Errors:", errors)

