# DNA project

# Imports
import random
import sys
import copy
import numpy as np
from direct_solvers import get_optimal_row_solution

# Globals:
DNA_BASES = 'ACGT'
DEBUG = False

# Functions to set up experiments
def create_random_strand(length):
    """
    Generates a random DNA strand of the given length.

    Parameters:
    length (int): The length of the DNA strand.

    Returns:
    str: A string representing the DNA strand consisting of 'A', 'C', 'G', 'T'.
    """
    return ''.join(random.choice("ACGT") for _ in range(length))


def create_dna_list(num_of_strands, length_of_strands):
    """
    Generates a list of DNA strands, to stimulate a row in the synthesis machine.

    Parameters:
    num_of_strands (int): The number of strands to generate.
    length_of_strands (int): The length of the DNA strands.

    Returns:
    list: A list of strings representing the DNA strands consisting of 'A', 'C', 'G', 'T'.
    """
    return [create_random_strand(length_of_strands) for _ in range(num_of_strands)]


# Functions for consuming
def consume_base(index, dna_list):
    if 0 <= index <= len(dna_list):
        dna_list[index] = dna_list[index][1:]
        if dna_list[index] == "":
            dna_list.pop(index)
    else:
        print("index is out of range!")


def consume_logic_random(dna_list, base):
    assert dna_list

    # Filter the strands that begin with the given base
    filtered_strands = [temp_strand for temp_strand in dna_list if temp_strand.startswith(base)]

    if len(filtered_strands) == 0:
        #        print("Debug. No base available!") TODO: remove
        return None

    return random.choice(filtered_strands)


def consume_logic_lookahead_one_v1(dna_list, base):
    assert dna_list

    # Filter the strands that begin with the given base
    filtered_strands = [temp_strand for temp_strand in dna_list if temp_strand.startswith(base)]

    if len(filtered_strands) == 0:
        #        print("Debug. No base available!") TODO: remove
        return None

    next_base = get_next_base(base)

    # Look for a strand with the next base as a second character
    for strand in filtered_strands:
        if len(strand) > 1 and strand[1] == next_base:
            return strand

    # Did not find any strand with the next base as the next character, returning a random strand
    return random.choice(filtered_strands)


# Can be written with better complexity
def consume_logic_lookahead_one_v2(dna_list, base):
    assert dna_list

    # Filter the strands that begin with the given base
    filtered_strands = [temp_strand for temp_strand in dna_list if temp_strand.startswith(base)]

    if len(filtered_strands) == 0:
        #        print("Debug. No base available!") TODO: remove
        return None

    next_base = get_next_base(base)
    shifted_dna_bases = shift_string(DNA_BASES, next_base)

    is_base_needed = []
    for curr_base in shifted_dna_bases[:-1]: # Last char will have special treatment
        is_base_needed.append(not any(s[0] == curr_base for s in dna_list))
    temp_counter = 0
    for strand in dna_list:
        if strand[0] == shifted_dna_bases[-1]:
            temp_counter += 1
            if temp_counter == 2:
                break
    is_base_needed.append(temp_counter > 1)

    for i, curr_base in enumerate(shifted_dna_bases):
        for strand in filtered_strands:
            if len(strand) > 1 and strand[1] == curr_base and (is_base_needed[i] is True):
                return strand

    # Did not find any strand with the next base as the next character, returning a random strand
    return random.choice(filtered_strands)

def consume_logic_lookahead_one_v3(dna_list, base):
    assert dna_list

    # Filter the strands that begin with the given base
    filtered_strands = [temp_strand for temp_strand in dna_list if temp_strand.startswith(base)]

    if len(filtered_strands) == 0:
        #        print("Debug. No base available!") TODO: remove
        return None

    bases_counter_d = {}
    for curr_base in DNA_BASES:
        bases_counter_d[curr_base] = 0

    for strand in dna_list:
        bases_counter_d[strand[0]] += 1

    # Removing 1 from current base, since it will be deleted by one anyway
    bases_counter_d[base] -= 1

    # Edge case: last base in last strand
    max_key = max(bases_counter_d, key=bases_counter_d.get)
    if max_key == 0:
        print("edge case!")
        return random.choice(filtered_strands)

    # Pick a strand that will unlock a base that is least common in current synthesizeable (funny word) indices
    sorted_by_value = sorted(bases_counter_d.items(), key=lambda x: x[1])
    for curr_base, _ in sorted_by_value:
        for strand in filtered_strands:
            if len(strand) > 1 and strand[1] == curr_base:
                return strand

    # Did not find any strand with the next base as the next character, returning a random strand
    return random.choice(filtered_strands)


def consume_with_logic(dna_list, base, func):

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
    index = DNA_BASES.find(current_base)

    if index == -1:
        raise ValueError(f"Base '{current_base}' is not a valid DNA base.")

    # Return the next base in the string, wrapping around cyclically
    next_index = (index + 1) % len(DNA_BASES)
    return DNA_BASES[next_index]


def shift_string(s, char):
    """Shifts string `s` to begin with the first occurrence of `char`."""
    index = s.find(char)

    if index == -1:
        raise ValueError(f"Character '{char}' not found in string.")

    return s[index:] + s[:index]  # Rotate the string


# DNA synthesiser machine TODO: currently only a single row is supported
def synthezise(dna_list, func, should_visualize = False):
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
    cycle_results = []
    percent_results = []
    for _ in range(100):
        cycles, percent = run_single_experiment(num_of_strands, len_of_strands, logic_func)
        cycle_results.append(cycles)
        percent_results.append(percent)

    mean = np.mean(cycle_results)
    percent_mean = np.mean(percent_results)
    print(f"The mean is: {mean}")
    print(f"The percent is is: {percent_mean}")


def run_single_experiment(num_of_strands, len_of_strands, logic_func, own_list=None):
    if own_list is None:
        dna_list = create_dna_list(num_of_strands, len_of_strands)
    else:
        dna_list = own_list

    dna_list_copy = copy.deepcopy(dna_list)
    dna_list_copy_2 = copy.deepcopy(dna_list)

    optimal_soultion = get_optimal_row_solution(dna_list_copy_2)
    print(f"The optimal solution is: {optimal_soultion}")


    num_of_cycles = synthezise(dna_list, logic_func)
    # print(f"The num of cycles for selected logic: {num_of_cycles}")

    num_of_cycles_control_group = synthezise(dna_list_copy, consume_logic_random)
    # print(f"The num of cycles for random logic: {num_of_cycles_control_group}")

    # Stats:
    diff_between_optimal_cycles_tested_logic = num_of_cycles - optimal_soultion
    diff_between_optimal_cycles_random_logic = num_of_cycles_control_group - optimal_soultion

    cycles_saved = num_of_cycles_control_group - num_of_cycles
    percentage_save = num_of_cycles_control_group / num_of_cycles



    print(f"The amount of cycles saved: {cycles_saved}")
    print(f"Diff from optimal solution: {diff_between_optimal_cycles_tested_logic}")

    return cycles_saved, percentage_save

import os
import sys
import time

def clear_console():
    # Works for both Windows and Linux/Mac
    os.system('cls' if os.name == 'nt' else 'clear')

if __name__ == '__main__':
    #run_experiment(4, 200, consume_logic_lookahead_one_v3)

    run_single_experiment(5, 5, consume_logic_lookahead_one_v3)

