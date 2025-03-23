import random
import math
import copy
import numpy as np
import matplotlib.pyplot as plt

# Ensure these are correctly exported in main.py
from main import create_random_strand, DNA_BASES
from main import consume_logic_lookahead_one_v1, consume_logic_lookahead_one_v2
from main import consume_logic_lookahead_one_v3, consume_logic_lookahead_one_v4, consume_logic_random
from main import consume_with_logic


###############################################################################
# Helper Functions: Histograms and Imbalance
###############################################################################

def compute_histogram(strand):
    """
    Compute the histogram of DNA bases for a given strand.

    :param strand: A string representing the DNA strand (e.g., "ACGTAC").
    :return: A tuple (count_A, count_C, count_G, count_T).
    """
    return (strand.count('A'),
            strand.count('C'),
            strand.count('G'),
            strand.count('T'))


def add_hist(h1, h2):
    """
    Perform element-wise addition of two 4-element histograms.

    :param h1: A 4-tuple of integers.
    :param h2: Another 4-tuple of integers.
    :return: A 4-tuple representing the element-wise sum.
    """
    return (h1[0] + h2[0],
            h1[1] + h2[1],
            h1[2] + h2[2],
            h1[3] + h2[3])


def imbalance_std(hist):
    """
    Compute the standard deviation of a 4-element histogram.

    :param hist: A 4-tuple (count_A, count_C, count_G, count_T).
    :return: A float indicating the imbalance (lower is more balanced).
    """
    total = sum(hist)
    if total == 0:
        return 0.0
    mean = total / 4.0
    return math.sqrt(((hist[0] - mean) ** 2 +
                      (hist[1] - mean) ** 2 +
                      (hist[2] - mean) ** 2 +
                      (hist[3] - mean) ** 2) / 4.0)


###############################################################################
# Balanced Histogram Allocation
###############################################################################

def assign_strands_balanced_aggregate(dna_strands, m):
    """
    Allocate DNA strands to m rows so that each row's overall base histogram
    is as balanced as possible.

    The algorithm works as follows:
      1. For each strand, compute its histogram and peak ratio (max(hist)/length).
      2. Sort strands by descending peak ratio (most skewed first).
      3. For each strand (in sorted order), assign it to the row (not yet full) that
         minimizes the increase in imbalance (standard deviation of cumulative histogram).

    :param dna_strands: List of DNA strands (total must be m*n).
    :param m: Number of rows.
    :return: A list of m rows, each row being a list of strands.
    """
    total = len(dna_strands)
    n = total // m
    strand_info = []

    for s in dna_strands:
        hist_s = compute_histogram(s)
        length_s = len(s)
        peak_s = max(hist_s)
        pr = peak_s / length_s if length_s > 0 else 0
        strand_info.append((s, hist_s, pr))

    # Sort by peak ratio descending
    strand_info.sort(key=lambda x: x[2], reverse=True)

    # Initialize rows and cumulative histograms
    Rows = [[] for _ in range(m)]
    RowHist = [(0, 0, 0, 0) for _ in range(m)]

    # Greedy assignment of each strand to the row with the smallest imbalance increase
    for (s, hist_s, pr) in strand_info:
        best_row = None
        best_imbalance = float('inf')
        for r in range(m):
            if len(Rows[r]) < n:
                candidate_hist = add_hist(RowHist[r], hist_s)
                candidate_imb = imbalance_std(candidate_hist)
                if candidate_imb < best_imbalance:
                    best_imbalance = candidate_imb
                    best_row = r
        Rows[best_row].append(s)
        RowHist[best_row] = add_hist(RowHist[best_row], hist_s)

    return Rows


###############################################################################
# Random Assignment (Baseline)
###############################################################################

def assign_strands_random(dna_strands, m):
    """
    Randomly assign DNA strands to m rows.

    :param dna_strands: List of DNA strands (total must be m*n).
    :param m: Number of rows.
    :return: A list of m rows, each containing n strands.
    """
    total = len(dna_strands)
    n = total // m
    idxs = list(range(total))
    random.shuffle(idxs)
    Rows = []
    start = 0
    for i in range(m):
        subset = idxs[start:start + n]
        row_strands = [dna_strands[j] for j in subset]
        Rows.append(row_strands)
        start += n
    return Rows


###############################################################################
# Synthesis Simulation (No Cycle Limit)
###############################################################################

def synthezise_no_limit(dna_list, consume_func):
    """
    Simulate synthesis of a single row without imposing a cycle limit.

    The synthesis machine cycles through DNA_BASES until all strands are synthesized.
    The consumption function determines which strand is advanced at each cycle.

    :param dna_list: List of DNA strands.
    :param consume_func: Function to choose which strand to consume (e.g., lookahead or random).
    :return: Total number of cycles required.
    """
    cycles = 0
    base_idx = 0
    while dna_list:
        base = DNA_BASES[base_idx]
        consume_with_logic(dna_list, base, consume_func)
        cycles += 1
        base_idx = (base_idx + 1) % len(DNA_BASES)
    return cycles


def simulate_matrix_synthesis(matrix, consume_func):
    """
    Simulate synthesis of an entire matrix of DNA strands with no cycle limit.

    Each row is synthesized independently; the total synthesis time is the maximum
    number of cycles among all rows.

    :param matrix: List of rows (each row is a list of DNA strands).
    :param consume_func: Consumption function used for synthesis.
    :return: Maximum number of cycles among the rows.
    """
    row_cycles = []
    for row in matrix:
        row_copy = copy.deepcopy(row)
        c = synthezise_no_limit(row_copy, consume_func)
        row_cycles.append(c)
    return max(row_cycles)


###############################################################################
# Experiment: Compare Allocation Methods
###############################################################################

def run_experiment(m, n, strand_length, num_trials=5):
    """
    Run multiple trials for a given matrix dimension (m rows, n columns).

    For each trial:
      - Generate m*n random DNA strands (each of length 'strand_length').
      - Allocate strands using both Balanced Histogram Allocation and Random Assignment.
      - Simulate synthesis using various consumption functions.
      - Compute the average cycle counts and the percentage improvement of Balanced vs Random.

    :param m: Number of rows.
    :param n: Number of strands per row.
    :param strand_length: Length of each DNA strand.
    :param num_trials: Number of trials for averaging.
    :return: Tuple (avg_results, improvement) where:
             avg_results: dict with average cycles for "Balanced" and "Random" methods.
             improvement: dict with percentage improvement per consumption function.
    """
    consume_funcs = {
        "Lookahead_v1": consume_logic_lookahead_one_v1,
        "Lookahead_v2": consume_logic_lookahead_one_v2,
        "Lookahead_v3": consume_logic_lookahead_one_v3,
        "Lookahead_v4": consume_logic_lookahead_one_v4,
        "Random": consume_logic_random
    }
    total_strands = m * n
    results = {
        "Balanced": {name: [] for name in consume_funcs},
        "Random": {name: [] for name in consume_funcs}
    }

    for trial in range(1, num_trials + 1):
        print(f"Trial {trial}/{num_trials} (m={m}, n={n})")
        strands = [create_random_strand(strand_length) for _ in range(total_strands)]
        balanced_rows = assign_strands_balanced_aggregate(strands, m)
        random_rows = assign_strands_random(strands, m)

        for name, func in consume_funcs.items():
            cycles_bal = simulate_matrix_synthesis(balanced_rows, func)
            cycles_rand = simulate_matrix_synthesis(random_rows, func)
            results["Balanced"][name].append(cycles_bal)
            results["Random"][name].append(cycles_rand)
            print(f"  {name}: Balanced = {cycles_bal} cycles, Random = {cycles_rand} cycles")
        print("-" * 40)

    avg_results = {"Balanced": {}, "Random": {}}
    improvement = {}
    for method in results:
        for name in results[method]:
            avg = np.mean(results[method][name])
            avg_results[method][name] = avg
    for name in consume_funcs:
        rand_avg = avg_results["Random"][name]
        bal_avg = avg_results["Balanced"][name]
        improvement[name] = ((rand_avg - bal_avg) / rand_avg) * 100 if rand_avg != 0 else 0.0

    print("\nAverage Synthesis Cycles:")
    for method in avg_results:
        print(f"  {method}:")
        for name, avg in avg_results[method].items():
            print(f"    {name}: {avg:.2f} cycles")
    print("\nPercentage Improvement of Balanced vs Random:")
    for name, imp in improvement.items():
        print(f"  {name}: {imp:.2f}%")

    return avg_results, improvement


###############################################################################
# Experiments Over Multiple Matrix Sizes with Enhanced Graphs
###############################################################################

def run_experiments_for_sizes(sizes, strand_length, num_trials):
    """
    For each square matrix size (m=n=size), run experiments and record:
      1) The percentage improvement (excluding the Random baseline).
      2) The absolute synthesis cycles for a chosen consumption function.
      3) A separate graph showing the difference (Random - Balanced).

    :param sizes: List of matrix sizes (integers).
    :param strand_length: Length of each DNA strand.
    :param num_trials: Number of trials per matrix size.
    :return: Tuple (improvement_data, avg_cycles_data_bal, avg_cycles_data_rand).
    """
    consume_names = ["Lookahead_v1", "Lookahead_v2", "Lookahead_v3", "Lookahead_v4", "Random"]
    consume_names_filtered = [name for name in consume_names if name != "Random"]

    improvement_data = {name: [] for name in consume_names_filtered}
    avg_cycles_data_bal = {name: [] for name in consume_names}
    avg_cycles_data_rand = {name: [] for name in consume_names}

    for size in sizes:
        m = size
        n = size
        print(f"\nRunning experiment for matrix size: {m} x {n}, strand_length={strand_length}")
        avg_results, improvement = run_experiment(m, n, strand_length, num_trials)

        for name in consume_names:
            avg_cycles_data_bal[name].append(avg_results["Balanced"][name])
            avg_cycles_data_rand[name].append(avg_results["Random"][name])
        for name in consume_names_filtered:
            improvement_data[name].append(improvement[name])

    # Plot (1): Percentage improvement vs. matrix size (excluding Random)
    plt.figure(figsize=(10, 6))
    for name in consume_names_filtered:
        plt.plot(sizes, improvement_data[name], marker='o', label=name)
    plt.xlabel("Matrix Size (m = n)")
    plt.ylabel("Percentage Improvement (%)")
    plt.title("Percentage Improvement of Balanced vs Random (No Cycle Limit)")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot (2): Absolute cycles for a chosen consumption function
    chosen = "Lookahead_v4"
    random_line = avg_cycles_data_rand[chosen]
    balanced_line = avg_cycles_data_bal[chosen]

    plt.figure(figsize=(10, 6))
    plt.plot(sizes, random_line, marker='s', linestyle='--', color='blue', label="Random")
    plt.plot(sizes, balanced_line, marker='o', linestyle='-', color='orange', label="Balanced")
    plt.xlabel("Matrix Size (m = n)")
    plt.ylabel("Average Synthesis Cycles")
    plt.title(f"Absolute Cycles vs Matrix Size for {chosen} (No Cycle Limit)")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot (3): Separate difference graph: (Random - Balanced)
    difference_line = [r - b for r, b in zip(random_line, balanced_line)]
    plt.figure(figsize=(10, 6))
    plt.plot(sizes, difference_line, marker='^', color='red')
    plt.xlabel("Matrix Size (m = n)")
    plt.ylabel("Difference in Cycles (Random - Balanced)")
    plt.title(f"Difference in Cycles for {chosen} (No Cycle Limit)")
    plt.grid(True)
    plt.show()

    return improvement_data, avg_cycles_data_bal, avg_cycles_data_rand


###############################################################################
# Main Runner
###############################################################################

if __name__ == '__main__':
    sizes = [4, 6, 8, 10, 12, 15, 20, 30]  # Square matrix sizes to test
    strand_length = 200
    num_trials = 10  # Trials per matrix size

    print("Running Balanced vs Random (No Cycle Limit) across multiple matrix sizes...")
    run_experiments_for_sizes(sizes, strand_length, num_trials)
