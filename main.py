# DNA project

# Imports
import random
import copy
import numpy as np
from graphs import plot_simple_graph
from lookahead import consume_logic_lookahead_one_v1, consume_logic_lookahead_one_v2
from lookahead import consume_logic_lookahead_one_v3, consume_logic_lookahead_one_v4
from lookahead import consume_logic_random

# Globals:
DNA_BASES = 'ACGT'
DEBUG = False

# Functions to set up experiments
def create_random_strand(length):
    """
    Generates a random DNA strand of the given length.

    Parameters:
    :param length:  The length of the DNA strand.
    :return: A string describing the DNA strand.
    """
    return ''.join(random.choice("ACGT") for _ in range(length))


def create_dna_list(num_of_strands, length_of_strands):
    """
    Generates a list of DNA strands, to stimulate a row in the synthesis machine.

    :param num_of_strands: The number of strands.
    :param length_of_strands: The length of the strands.
    :return: A list of strings representing the DNA strands consisting of 'A', 'C', 'G', 'T'.
    """
    return [create_random_strand(length_of_strands) for _ in range(num_of_strands)]


# Functions for consuming
def consume_base(index, dna_list):
    """
    Remove a base from a dna strand, dna_list[index]

    :param index: The strand to pick from dna_list.
    :param dna_list: The list of the dna strands.
    """
    if 0 <= index <= len(dna_list):
        dna_list[index] = dna_list[index][1:]
        if dna_list[index] == "":
            dna_list.pop(index)
    else:
        print("index is out of range!")


def consume_with_logic(dna_list, base, func):
    """
    Use a consume function to decide which base to remove

    :param dna_list: The list of the dna strands.
    :param base: The base to consume
    :param func: The logic function to consume with.
    """

    assert dna_list

    selected_strand = func(dna_list, base)

    if selected_strand is None:
        return

    # Find the index of the selected strand in the original list
    selected_index = dna_list.index(selected_strand) if selected_strand else None

    consume_base(selected_index, dna_list)

    return


# Functions for DNA arithmetics

def get_next_base(current_base):
    """
    Remove a base from a dna strand, dna_list[index]

    :param current_base: The current base.
    :return: The next machine base, based on DNA_BASES.
    """
    index = DNA_BASES.find(current_base)

    if index == -1:
        raise ValueError(f"Base '{current_base}' is not a valid DNA base.")

    # Return the next base in the string, wrapping around cyclically
    next_index = (index + 1) % len(DNA_BASES)
    return DNA_BASES[next_index]


def shift_string(s, char):
    """
    Shift a string to start with the char given. Used in lookahead for efficiency.
    :param s: The string to shift
    :param char: The char that will be the start of the string
    :return: The shifted string.
    """
    index = s.find(char)

    if index == -1:
        raise ValueError(f"Character '{char}' not found in string.")

    return s[index:] + s[:index]  # Rotate the string


def synthesize(dna_list, func, should_visualize = False):
    """
    Given a dna_list that represents a row, synthesize the strands
    :param dna_list: The dna strands in a row
    :param func: The logic func that will decide at each step which strand to pick
    :param should_visualize: If true, and run through command line (and not some IDE)
    will interactively show the process
    :return: The amount of cycles it took to synthesize
    """
    base = DNA_BASES[0]
    num_of_cycles = 0
    while dna_list:
        if should_visualize:
            clear_console()
            for strand in dna_list:
                print(strand)
            time.sleep(0.1)
        consume_with_logic(dna_list, base, func)

        num_of_cycles += 1
        base = get_next_base(base)

    return num_of_cycles


def run_experiment(num_of_strands, len_of_strands, logic_func):
    """
    Run an experiment. Will generate a list of strands given the parameters, and synthesize them.
    One time it will synthesize with the random logic, the second time with the logic func given.
    :param num_of_strands: The amount of strands.
    :param len_of_strands: The len of each strand.
    :param logic_func: The logic function to test against the random picking algorithm.
    :return: The mean and percent of cycles saved (since this is repeated 300 times, we return a mean)
    """
    cycle_results = []
    percent_results = []
    for _ in range(300):
        cycles, percent = run_single_experiment(num_of_strands, len_of_strands, logic_func)
        cycle_results.append(cycles)
        percent_results.append(percent)

    mean = np.mean(cycle_results)
    percent_mean = np.mean(percent_results)
    print(f"The mean is: {mean}")
    print(f"The percent is is: {percent_mean}")

    return mean, percent_mean


def run_single_experiment(num_of_strands, len_of_strands, logic_func, own_list=None):
    """
    An auxiliary function of run_experiment. Runs a single experiment.
    :param num_of_strands: How many strands in the row.
    :param len_of_strands: The len of each strand.
    :param logic_func: The function to use to consume.
    :param own_list: A list of DNA strands. If given, the first two parameters are ignored, and the list will be used.
    :return: The difference in cycles between random and logic_func (also return the percent change).
    """
    if own_list is None:
        dna_list = create_dna_list(num_of_strands, len_of_strands)
    else:
        dna_list = own_list

    dna_list_copy = copy.deepcopy(dna_list)

    num_of_cycles = synthesize(dna_list, logic_func, should_visualize=False)
    # print(f"The num of cycles for selected logic: {num_of_cycles}")

    num_of_cycles_control_group = synthesize(dna_list_copy, consume_logic_random, should_visualize=False)
    # print(f"The num of cycles for random logic: {num_of_cycles_control_group}")

    cycles_saved = num_of_cycles_control_group - num_of_cycles
    percentage_save = num_of_cycles_control_group / num_of_cycles

    return cycles_saved, percentage_save


logic_functions = [
        consume_logic_lookahead_one_v1,
        consume_logic_lookahead_one_v2,
        consume_logic_lookahead_one_v3,
        consume_logic_lookahead_one_v4]

# Get function names as strings
logic_function_names = [func.__name__ for func in logic_functions]


def test_all_heuristics(num_of_strands_list):
    """
    Tests all of our heuristics at the same time.
    :param num_of_strands_list: Each cell in the list will run an experiment with the cell's number of strands.
    """
    import matplotlib.pyplot as plt

    num_of_logic_funcs = len(logic_functions)

    cycles_lists_mean = [[] for _ in range(num_of_logic_funcs)]
    cycles_lists_percent = [[] for _ in range(num_of_logic_funcs)]

    for j in range(len(num_of_strands_list)):
        cycles_lists = [[] for _ in range(num_of_logic_funcs)]
        # Repeat experiment 100 times
        print(f'starting cycle number {j}')
        cycles_saved_lists = [[] for _ in range(num_of_logic_funcs)]
        percent_saved_lists = [[] for _ in range(num_of_logic_funcs)]
        for _ in range (50):
            dna_list = create_dna_list(num_of_strands_list[j], 50)
            for i in range(num_of_logic_funcs):
                dna_list_copy = copy.deepcopy(dna_list)
                cycles_saved, percent_saved = run_single_experiment(_,_, logic_functions[i], dna_list_copy)
                cycles_saved_lists[i].append(cycles_saved)
                percent_saved_lists[i].append(percent_saved)
        for i in range(num_of_logic_funcs):
            cycles_lists_mean[i].append(np.mean(cycles_saved_lists[i]))
            cycles_lists_percent[i].append(np.mean(percent_saved_lists[i]))

    # Plot the graph
    plt.figure(figsize=(8, 5))
    color_list = ['b', 'g', 'r', 'm']
    for i in range(len(logic_function_names)):
        plt.plot(num_of_strands_list, cycles_lists_mean[i], label=logic_function_names[i], color=color_list[i])

    # Labels and title
    plt.xlabel("Number of Strands")
    plt.ylabel("Reduced Overhead Cycles")
    plt.title("Impact of Heuristic Choice on Overhead Cycle Reduction")
    plt.legend()  # Show legend
    plt.grid(True)  # Enable grid

    # Show the plot
    plt.show()

    #percent graph
    plt.figure(figsize=(8, 5))
    for i in range(len(logic_function_names)):
        plt.plot(num_of_strands_list, cycles_lists_percent[i], label=logic_function_names[i], color=color_list[i])

    # Labels and title
    plt.xlabel("Number of Strands")
    plt.ylabel("Reduced Overhead Cycles (%)")
    plt.title("Impact of Heuristic Choice on Overhead Cycle Reduction (%)")
    plt.legend()  # Show legend
    plt.grid(True)  # Enable grid

    plt.show()

import os
import sys
import time

def clear_console():
    # Works for both Windows and Linux/Mac
    os.system('cls' if os.name == 'nt' else 'clear')


def shift_dna_bases(strand: str, offset: int) -> str:
    """
    Shifts each base in the DNA strand by the given offset in the order A -> C -> G -> T -> A.

    :param strand: A string representing the DNA sequence (e.g., "ACGT").
    :param offset: An integer specifying how many steps to shift each base.
    :return: A new DNA strand with bases shifted.
    """
    base_map = {base: i for i, base in enumerate(DNA_BASES)}

    shifted_strand = "".join(DNA_BASES[(base_map[base] + offset) % 4] for base in strand)
    return shifted_strand


def count_dna_bases(strand: str) -> dict:
    """
    Counts the number of occurrences of each DNA base (A, C, G, T) in the given strand.

    :param strand: A string representing the DNA sequence (e.g., "ACGT").
    :return: A dictionary with counts of each base.
    """
    return {base: strand.count(base) for base in DNA_BASES}


def balance_dna_strands(strand1: str, strand2: str) -> tuple:
    """
    Determines the optimal shift for each strand to achieve the best balance of bases.

    :param strand1: The first DNA strand.
    :param strand2: The second DNA strand.
    :return: A tuple containing the best shift values for strand1 and strand2.
    """
    best_shift1, best_shift2 = 0, 0
    min_imbalance = float('inf')
    dna_strand_length = len(strand1)

    for shift1 in range(len(DNA_BASES)):
        shifted1 = shift_dna_bases(strand1, shift1)
        count1 = count_dna_bases(shifted1)

        for shift2 in range(len(DNA_BASES)):
            shifted2 = shift_dna_bases(strand2, shift2)
            count2 = count_dna_bases(shifted2)

            total_balance = sum(abs(count1[base] + count2[base] - (2*(dna_strand_length / 4))) for base in DNA_BASES)

            if total_balance < min_imbalance:
                min_imbalance = total_balance
                best_shift1, best_shift2 = shift1, shift2

    best_shifted_strand1 = shift_dna_bases(strand1, best_shift1) + DNA_BASES[best_shift1]
    best_shifted_strand2 = shift_dna_bases(strand2, best_shift2) + DNA_BASES[best_shift2]

    # best_shifted_strand1 = shift_dna_bases(strand1, best_shift1)
    # best_shifted_strand2 = shift_dna_bases(strand2, best_shift2)

    return best_shifted_strand1, best_shifted_strand2


def balance_dna_strands2(dna_strands: list) -> list:
    """
    Determines the optimal shift for each DNA strand to achieve the best balance of bases across all strands.

    :param dna_strands: Multiple DNA strands.
    :return: A tuple containing the best shifted strands.
    """
    import itertools

    best_shifts = [0] * len(dna_strands)
    min_imbalance = float('inf')
    dna_strand_length = len(dna_strands[0])

    # Iterate over all strands
    for shifts in itertools.product(range(len(DNA_BASES)), repeat=len(dna_strands)):
        total_balance = 0

        # For each strand, apply the shifts and compute the imbalance
        shifted_strands = []
        for idx, shift in enumerate(shifts):
            shifted_strand = shift_dna_bases(dna_strands[idx], shift)
            count = count_dna_bases(shifted_strand)
            shifted_strands.append(shifted_strand)

            # Add up imbalance for this strand
            total_balance += sum(abs(count[base] - (dna_strand_length / 4)) for base in DNA_BASES)


        # Update the best shifts if this combination reduces imbalance
        if total_balance < min_imbalance:
            min_imbalance = total_balance
            best_shifts = list(shifts)

    # Apply the best shifts to each strand and return them
    best_shifted_strands = [shift_dna_bases(strand, best_shifts[idx]) + DNA_BASES[best_shifts[idx]] for idx, strand in
                            enumerate(dna_strands)]

    return best_shifted_strands

def process_dna_pairs(dna_list: list) -> list:
    """
    Processes pairs of DNA strands from the list, applying the optimal shift and appending
    the corresponding base to each strand.

    :param dna_list: A list of DNA strands (even size).
    :return: A list of processed DNA strands with shift indicators appended.
    """
    if len(dna_list) % 2 != 0:
        raise ValueError("The list must have an even number of DNA strands.")

    processed_strands = []

    for i in range(0, len(dna_list), 2):
        strand1 = dna_list[i]
        strand2 = dna_list[i + 1]

        best_shifted_strand1, best_shifted_strand2 = balance_dna_strands(strand1, strand2)

        # Append both shifted strands to the result list
        processed_strands.append(best_shifted_strand1)
        processed_strands.append(best_shifted_strand2)

    return processed_strands


def run_experiment_change_dna_strands(num_of_strands, len_of_strands, logic_func):
    """
    Runs a single experiment of that involves changing DNA strands.
    :param num_of_strands: The number of strands to use.
    :param len_of_strands: The length of each strand
    :param logic_func: Which logic func to use to consume
    :return: The mean and percent of change with respect to random.
    """
    cycle_results = []
    percent_results = []
    for _ in range(300):
        dna_list = create_dna_list(num_of_strands, len_of_strands)
        shifted_list = process_dna_pairs(dna_list)

        num_of_cycles = synthesize(shifted_list, logic_func, should_visualize=False)
        # print(f"The num of cycles for selected logic: {num_of_cycles}")

        num_of_cycles_control_group = synthesize(dna_list, logic_func, should_visualize=False)
        # print(f"The num of cycles for random logic: {num_of_cycles_control_group}")

        cycles_saved = num_of_cycles_control_group - num_of_cycles
        percentage_save = num_of_cycles_control_group / num_of_cycles

        cycle_results.append(cycles_saved)
        percent_results.append(percentage_save)

    mean = np.mean(cycle_results)
    percent_mean = np.mean(percent_results)

    return mean, percent_mean

if __name__ == '__main__':
    test_all_heuristics([2,3,4,5,6,7,8,10,12,15,20,25,30])
