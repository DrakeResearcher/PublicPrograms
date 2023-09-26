#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import re

from decimal import *
from mpl_toolkits.axes_grid1.inset_locator import  mark_inset

#TODO: isolate x=24 to 30 and add it in an inset

def main():
    diffs = calculate_differences(
        [ # this list contains all files to be read
            'data/truncated/111SPOW.MAT',
            'data/not-truncated/111SPOW.MAT'
        ], 
        'not-truncated' # baseline E value
    )

    chart_values(diffs, True)

# Does the actual charting
def chart_values(values: [[int, Decimal]], show=False):
    colors = [
        "black",
        "red",
        "blue",
        "orange",
        "purple",
        "green",
        "yellow"
    ]

    styles = [
        "-s",
        "-^",
        "-o",
        "-x",
        "-v",
        "-S",
    ]

    fig, main = plt.subplots() # get figure and axis

    for index, data_set in enumerate(values): # iterate thru values with index
        raw_data = list(map(list, zip(*data_set)))  # transpose data_set:
                                                    # goes from [[num, E-diff]_1, []_2, ..., []_n]
                                                    # to [[num_1, num_2, ..., num_n], [Ediff_1, ..., Ediff_n]]
                                                    # (all of our x values and all of our y values)

        main.plot(
            raw_data[0], raw_data[1], # x values, y values
            styles[index], # solid line with circle markers
            markersize=2,
            color=colors[index],
        )

    dims = [0.45, 0.45, 0.4, 0.4] # dimensions of inset graph: x, y, w, h
    sub = fig.add_axes(dims) # make new axes with dimensions

    for index, data_set in enumerate(values): # same as above except we will crop this graph later
        raw_data = list(map(list, zip(*data_set)))

        sub.plot(
            raw_data[0], raw_data[1],
            styles[index],
            markersize=2,
            color=colors[index],
        )
   
    main.legend(["Truncated", "Non-Truncated"], frameon=False, loc="lower left")
    main.set_xlabel("Number of Terms")
    main.set_ylabel("Fractional Energy Difference")
    main.set_yscale("log")
    main.set_title("Convergence Comparison")

    main.tick_params(
        top=True,
        right=True,
        direction="in"
    )

    sub.set_yscale("log")
    sub.set_xlim([6000, 9_500]) # crop the inset graph
    sub.set_ylim([1e-25, 5e-23])

    # sub.tick_params(
    #     left=False, labelleft=False,
    #     bottom=False, labelbottom=False
    # )

    # sub.yaxis.set_minor_locator(ticker.NullLocator())

    mark_inset(main, sub, loc1=1, loc2=3, fc="none", ec='0.5')

    # plt.savefig("graph.pgf", dpi=300)

    if show:
        plt.show()


def calculate_differences(filenames: [str], dir: str):
    e_values = []

    for filename in filenames:
        e_values.append(compile_e_values(filename))

    file = open(f'data/{dir}/111SPOW.MAT') # file for E_extp
    e_extp_line = file.readlines().pop() # get final line
    file.close()

    e_extp = Decimal(generate_values(e_extp_line)[1])

    result = [[] for _ in range(len(filenames))] # initialize result to contain n lists where n = number of files

    for index, values in enumerate(e_values): # loop through all energy value lists with index
        for value in values: # for each energy value in the current list
            result[index].append([
                value[0], # number of terms
                (e_extp - value[1]) / e_extp 
            ])

    return result


# Does the bulk of the work. Reads main file and compiles data
# Return value: [ [# terms, E]1, []2, ..., []n ]
# For example, 111S251's entry would be: [ 251, abs(E_111S251 - E_extp) / E_extp ]
def compile_e_values(filename: str):
    file = open(filename, "r")
    lines = file.readlines()
    file.close()


    e_values = []
    e_values_converted = []

    i = 1
    for line in lines: # locate below which the energy values lie
        if "ENERGIES" in line:
            break
        i += 1

    lines = lines[i:]
    lines.pop() # remove final line of file which contains value for e_extp

    for line in lines:
        e_values.append(generate_values(line))

    for value in e_values:
        e_values_converted.append([int(value[0]), Decimal(value[1])])

    return e_values_converted

# Extracts name and energy value from a given line.
# If search is true, also calls get_omega to get the relevant omega value
def generate_values(line: str):
    trimmed_and_split = line.strip().split()
    return [trimmed_and_split[0], trimmed_and_split[1]]

def get_omega(wave: str, filename: str):
    white_space_regex = r"(?<!\()[ ]+" # find all whitespace that isnt preceded by an opening parenthesis.
                                       # necessary because some lines have a value such as 13( 3) in which
                                       # the space within the parens cause the value to be erroneously split,
                                       # which is bad because the number of slices the line is split into is
                                       # how we determine whether it is doubled or tripled. 

    file = open(filename, 'r')
    trunc = False

    if 'POW' in filename: # useless for now, will probably be useful later
        ext = 'POW'
    else:
        ext = 'POL'

    pow_file = f"111S{wave}.{ext}"

    for line in file.readlines():
        if pow_file in line:
            line_arr = re.split(white_space_regex, line.strip())
            
            if len(line_arr) >= 7: # should only be 6 or 7. <=6 -> doubled, >=7 -> tripled
                triple = True
            else:
                triple = False

            if int(line_arr[3]) > 0:
                trunc = True

    file.close()

    path = build_path(pow_file, triple, trunc)

    file = open(path, 'r')

    for i in range(3): 
        line = file.readline()

    line_arr = re.split(white_space_regex, line.strip())

    file.close()

    return line_arr[4]
    
def build_path(base: str, triple: bool, trunc: bool):
    path = '../wave/'

    if base[-3:] == 'POW':
        inf = True
    else:
        inf = False

    if inf:
        path += 'infinite-mass/'
    else:
        path += 'finite-mass/'

    if triple:
        path += 'tripled/'
    else:
        path += 'doubled/'
    
    if trunc:
        path += 'truncated/'
    else:
        path += 'not-truncated/'

    path += base

    return path

if (__name__ == "__main__"):
    main()