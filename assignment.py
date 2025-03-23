import random
import math
import copy
import numpy as np
import matplotlib.pyplot as plt

# Make sure these are correctly exported in your main.py
from main import create_random_strand, DNA_BASES
from main import consume_logic_lookahead_one_v1, consume_logic_lookahead_one_v2
from main import consume_logic_lookahead_one_v3, consume_logic_lookahead_one_v4, consume_logic_random
from main import consume_with_logic


###############################################################################
# 1. Helper Functions for Histograms and Imbalance
###############################################################################

def compute_histogram(strand):
    """
    Computes the histogram of DNA bases for a given strand.

    :param strand: A string representing the DNA strand (consisting of 'A','C','G','T').
    :return: A 4-element tuple (count_A, count_C, count_G, count_T).
    """
    return (strand.count('A'),
            strand.count('C'),
            strand.count('G'),
            strand.count('T'))


def add_hist(h1, h2):
    """
    Performs element-wise addition of two 4-element histograms.

    :param h1: A tuple of 4 integers (histogram).
    :param h2: Another tuple of 4 integers (histogram).
    :return: A tuple representing the sum of h1 and h2, element by element.
    """
    return (h1[0] + h2[0],
            h1[1] + h2[1],
            h1[2] + h2[2],
            h1[3] + h2[3])


def imbalance_std(hist):
    """
    Computes the standard deviation of the 4 histogram values.

    :param hist: A 4-element tuple representing (count_A, count_C, count_G, count_T).
    :return: A float indicating how imbalanced the histogram is (lower = more balanced).
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
# 2. Balanced Histogram Allocation Algorithm
###############################################################################

def assign_strands_balanced_aggregate(dna_strands, m):
    """
    Assigns strands to m rows using the Balanced Histogram Allocation method.

    This algorithm sorts all strands by their 'peak ratio' (max base-count / length),
    then assigns each strand to the row that, when the strand's histogram is added,
    yields the smallest increase in row imbalance (measured by standard deviation).

    :param dna_strands: A list of DNA strands to be allocated. Total must be m*n.
    :param m: Number of rows in the matrix.
    :return: A list (length = m) of row allocations (each row is a list of strands).
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

    # Sort strands by peak ratio descending (more skewed first)
    strand_info.sort(key=lambda x: x[2], reverse=True)

    # Initialize m empty rows and cumulative histograms
    Rows = [[] for _ in range(m)]
    RowHist = [(0, 0, 0, 0) for _ in range(m)]

    # Greedy assignment of strands
    for (s, hist_s, pr) in strand_info:
        best_row = None
        best_imbalance = float('inf')
        for r in range(m):
            if len(Rows[r]) < n:
                new_hist = add_hist(RowHist[r], hist_s)
                new_imb = imbalance_std(new_hist)
                if new_imb < best_imbalance:
                    best_imbalance = new_imb
                    best_row = r
        Rows[best_row].append(s)
        RowHist[best_row] = add_hist(RowHist[best_row], hist_s)

    return Rows


###############################################################################
# 3. Random Assignment (Baseline)
###############################################################################

def assign_strands_random(dna_strands, m):
    """
    Randomly assigns strands to m rows.

    :param dna_strands: A list of strands (total must be m*n).
    :param m: Number of rows in the matrix.
    :return: A list of m rows (lists of strands), each containing n strands.
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
# 4. Synthesis Simulation: No Cycle Limit
###############################################################################

def synthezise_no_limit(dna_list, consume_func):
    """
    Synthesizes a single row (dna_list) with no cycle limit.

    The machine cycles through DNA_BASES repeatedly until all strands are empty.
    The chosen consumption function is used to pick which strand to consume at each step.

    :param dna_list: A list of DNA strands.
    :param consume_func: A function that picks which strand to consume based on the current base.
    :return: The total number of cycles used once all strands are synthesized.
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
    Simulates synthesis on the entire matrix (list of rows) with no cycle limit.

    :param matrix: A list of rows (each row is a list of strands).
    :param consume_func: The consumption function used at each cycle (lookahead v1-v4 or random).
    :return: The maximum number of cycles among all rows (the matrix's total completion time).
    """
    row_cycles = []
    for row in matrix:
        row_copy = copy.deepcopy(row)
        c = synthezise_no_limit(row_copy, consume_func)
        row_cycles.append(c)
    return max(row_cycles)


###############################################################################
# 5. Experiment Function: Compare Balanced vs. Random Assignment
###############################################################################

def run_experiment(m, n, strand_length, num_trials=5):
    """
    Runs multiple trials for a single matrix dimension (m rows, n columns).
    Each trial:
      - Generates m*n random strands (each of length 'strand_length').
      - Allocates them via Balanced Histogram Allocation and via Random Assignment.
      - Synthesizes each matrix with each consumption function.
      - Computes average cycles over trials and the percentage improvement of Balanced over Random.

    :param m: Number of rows in the matrix.
    :param n: Number of strands per row.
    :param strand_length: Length of each DNA strand (in bases).
    :param num_trials: Number of trials to average.
    :return: (avg_results, improvement) where:
         avg_results is a dict of average cycles for "Balanced" and "Random",
         improvement is a dict with percentage improvement of Balanced vs Random.
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
        # Generate random strands
        strands = [create_random_strand(strand_length) for _ in range(total_strands)]

        # Balanced allocation
        balanced_rows = assign_strands_balanced_aggregate(strands, m)
        # Random allocation
        random_rows = assign_strands_random(strands, m)

        # Synthesize each matrix for every consumption function
        for name, func in consume_funcs.items():
            cycles_bal = simulate_matrix_synthesis(balanced_rows, func)
            cycles_rand = simulate_matrix_synthesis(random_rows, func)
            results["Balanced"][name].append(cycles_bal)
            results["Random"][name].append(cycles_rand)
            print(f"  {name}: Balanced = {cycles_bal} cycles, Random = {cycles_rand} cycles")
        print("-" * 40)

    # Compute average cycles
    avg_results = {"Balanced": {}, "Random": {}}
    improvement = {}

    for method in results:
        for name in results[method]:
            avg = np.mean(results[method][name])
            avg_results[method][name] = avg

    # Compute percentage improvement: ((Random - Balanced) / Random) * 100%
    for name in consume_funcs:
        rand_avg = avg_results["Random"][name]
        bal_avg = avg_results["Balanced"][name]
        if rand_avg != 0:
            improvement[name] = ((rand_avg - bal_avg) / rand_avg) * 100
        else:
            improvement[name] = 0.0

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
# 6. Testing Over Multiple Matrix Sizes with Enhanced Graphs
###############################################################################

def run_experiments_for_sizes(sizes, strand_length, num_trials):
    """
    For each matrix size (square: m=n=size),
    runs the experiment and records the percentage improvement for each consumption function.
    Plots graphs for percentage improvement and the difference in cycles for a chosen function.
    The "Random" consumption function is removed from the percentage improvement graph.

    :param sizes: List of integers for matrix sizes.
    :param strand_length: DNA strand length (in bases).
    :param num_trials: Number of trials per size.
    :return: Tuple (improvement_data, avg_cycles_data_bal, avg_cycles_data_rand).
    """
    consume_names = ["Lookahead_v1", "Lookahead_v2", "Lookahead_v3", "Lookahead_v4", "Random"]
    # For percentage improvement, remove "Random" since it's always 0.
    consume_names_filtered = [name for name in consume_names if name != "Random"]

    improvement_data = {name: [] for name in consume_names_filtered}
    avg_cycles_data_bal = {name: [] for name in consume_names}
    avg_cycles_data_rand = {name: [] for name in consume_names}

    for size in sizes:
        m = size
        n = size  # square matrix
        print(f"\nRunning experiment for matrix size: {m} x {n}, strand_length={strand_length}")
        avg_results, improvement = run_experiment(m, n, strand_length, num_trials)
        for name in consume_names:
            avg_cycles_data_bal[name].append(avg_results["Balanced"][name])
            avg_cycles_data_rand[name].append(avg_results["Random"][name])
        for name in consume_names_filtered:
            improvement_data[name].append(improvement[name])

    # Plot percentage improvement vs matrix size (excluding "Random")
    plt.figure(figsize=(10, 6))
    for name in consume_names_filtered:
        plt.plot(sizes, improvement_data[name], marker='o', label=name)
    plt.xlabel("Matrix Size (m = n)")
    plt.ylabel("Percentage Improvement (%)")
    plt.title("Percentage Improvement of Balanced vs Random (No Cycle Limit)")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Show difference in cycles for a chosen consumption function (e.g., Lookahead_v4)
    chosen = "Lookahead_v4"
    cycles_diff = [avg_cycles_data_rand[chosen][i] - avg_cycles_data_bal[chosen][i] for i in range(len(sizes))]
    plt.figure(figsize=(10, 6))
    plt.plot(sizes, cycles_diff, marker='o', color='red')
    plt.xlabel("Matrix Size (m = n)")
    plt.ylabel("Difference in Cycles (Random - Balanced)")
    plt.title(f"Difference in Synthesis Cycles for {chosen}")
    plt.grid(True)
    plt.show()

    # Plot absolute cycles for the chosen function
    plt.figure(figsize=(10, 6))
    plt.plot(sizes, avg_cycles_data_rand[chosen], marker='s', linestyle='--', label="Random")
    plt.plot(sizes, avg_cycles_data_bal[chosen], marker='o', linestyle='-', label="Balanced")
    plt.xlabel("Matrix Size (m = n)")
    plt.ylabel("Average Synthesis Cycles")
    plt.title(f"Absolute Cycles vs Matrix Size for {chosen}")
    plt.legend()
    plt.grid(True)
    plt.show()

    return improvement_data, avg_cycles_data_bal, avg_cycles_data_rand


###############################################################################
# 7. Main Runner
###############################################################################

if __name__ == '__main__':
    # Matrix sizes to test (square matrices)
    sizes = [4, 6, 8, 10, 12, 15, 20, 30]
    strand_length = 200
    num_trials = 10  # Number of trials per matrix size

    print("Running Balanced vs Random (No Cycle Limit) across multiple matrix sizes...")
    run_experiments_for_sizes(sizes, strand_length, num_trials)
