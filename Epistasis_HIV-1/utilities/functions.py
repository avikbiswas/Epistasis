from encodings.aliases import aliases
import re
import numpy as np
import pandas as pd
import math


# Read the reduced dictionary from a file and apply an offset to positions
def get_redu_dict(redux_file, offset):
    redux_dict = {}
    with open(redux_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            position = int(parts[0]) + offset
            groups = parts[1:]
            for i, group in enumerate(groups):
                key = (position, chr(65 + i))  # 'A', 'B', 'C', 'D', etc.
                redux_dict[key] = list(group.replace('-', ''))  # Remove '-' from the group
    return redux_dict


def unreduced_to_reduced(redux_dict, unreduced_pair):
    # Parse the unreduced pair
    wildtype, pos, mutate = unreduced_pair[0], int(unreduced_pair[1:-1]), unreduced_pair[-1]
    # Look up the reduced wildtype and mutate in the dictionary
    reduced_wildtype = next((key[1] for key, values in redux_dict.items() if key[0] == pos and wildtype in values), '-')
    reduced_mutate = next((key[1] for key, values in redux_dict.items() if key[0] == pos and mutate in values), '-')
    # Combine the reduced letters with the position
    reduced_pair = f"{reduced_wildtype}{pos}{reduced_mutate}"
    return reduced_pair


# Split a string of two pairs separated by a hyphen into a tuple
def split_pairs(pair_string):
    return tuple(pair.strip() for pair in pair_string.split('-') if pair.strip())


# Split a single pair like 'C140D' into its components
def split_pair(pair):
    wildtype = pair[0]
    position = int(pair[1:-1])  # Ensure position is an integer
    mutate = pair[-1]
    return wildtype, position, mutate


def DDE_dict(kn_file):
    dde_dict = {}
    with open(kn_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            columns = re.split(r'[,\t\s]+', line)
            pairs = columns[0]
            dde = float(columns[1])
            dde_dict[pairs] = dde
    return dde_dict
# get de from kn.her.all dictionary

#input: kn.her2.all
#output dde sorted list of tuples (pos1, pos2, ddE)
def kn_file(file_path):
    delta_delta_e_per_pair = []
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            columns = re.split(r'[,\t\s]+', line)
            mutations = columns[0]
            ddE = float(columns[1])
            positions = [int(pos) for pos in re.findall(r'\d+', mutations)]
            if len(positions) != 2:
                continue

            pos1, pos2 = positions
            delta_delta_e_per_pair.append((pos1, pos2, ddE))

    # Sort by ddE values
    delta_delta_e_per_pair_sorted = sorted(delta_delta_e_per_pair, key=lambda x: x[-1], reverse=True)
    return delta_delta_e_per_pair_sorted

def get_pos1_pos2(two_pairs):
    pair1 = split_pairs(two_pairs)[0]
    pair1_pos = get_pos(pair1)
    pair2 = split_pairs(two_pairs)[1]
    pair2_pos = get_pos(pair2)
    return int(pair1_pos),int(pair2_pos)



def unreduced_to_reduced_pair(redux_dict,unreduced_pair):
    wildtype1, pos1, mutate1 = split_pair(unreduced_pair.split('-')[0])
    wildtype2, pos2, mutate2 = split_pair(unreduced_pair.split('-')[1])

    reduced_wildtype1 = next((key[1] for key, values in redux_dict.items() if key[0] == pos1 and wildtype1 in values), '-')
    reduced_mutate1 = next((key[1] for key, values in redux_dict.items() if key[0] == pos1 and mutate1 in values), '-')

    reduced_wildtype2 = next((key[1] for key, values in redux_dict.items() if key[0] == pos2 and wildtype2 in values), '-')
    reduced_mutate2 = next((key[1] for key, values in redux_dict.items() if key[0] == pos2 and mutate2 in values), '-')

    reduced_pair = f"{reduced_wildtype1}{pos1}{reduced_mutate1}-{reduced_wildtype2}{pos2}{reduced_mutate2}"
    return reduced_pair


def get_wt(single_mutation):
    wt = single_mutation[:1]
    return wt


def get_pos(single_mutation):
    pos = int(single_mutation[1:-1])
    return pos

def read_list_from_file(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]
