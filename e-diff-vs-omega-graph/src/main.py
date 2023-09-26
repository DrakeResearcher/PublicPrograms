import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib
import re

from decimal import *
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)

#TODO: isolate x=24 to 30 and add it in an inset

def main():
    diffs = calculate_differences(
        [
            'data/truncated/111SPOW.MAT',
            'data/not-truncated/111SPOW.MAT'
        ], 
        'not-truncated'
    )

    chart_values(diffs, True)

# Does the actual charting
def chart_values(values: [[int, Decimal]], show=False):
    colors = [
        "red",
        "blue",
        "orange",
        "purple",
        "green",
        "yellow"
    ]

    fig, ax1 = plt.subplots()

    for index, data_set in enumerate(values):
        raw_data = list(map(list, zip(*data_set)))

        ax1.plot(
            raw_data[0], raw_data[1],
            "-",
            color=colors[index],
        )

    x, y, w, h = [0.40, 0.45, 0.4, 0.4]
    ax2 = fig.add_axes([x, y, w, h])

    for index, data_set in enumerate(values):
        raw_data = list(map(list, zip(*data_set)))

        ax2.plot(
            raw_data[0], raw_data[1],
            "-",
            color=colors[index],
        )
   
    ax1.legend(["Truncated", "Non-Truncated"], frameon=False, loc="lower left")
    ax1.set_xlabel("Number of Terms")
    ax1.set_ylabel("Fractional Energy Difference")
    ax1.set_yscale("log")
    ax1.set_title(r'$\frac{%.2f}{2}$' % (x))

    ax1.tick_params(
        top=True,
        right=True,
        direction="in"
    )

    ax2.set_yscale("log")
    ax2.set_xlim([6000, 9_200])
    ax2.set_ylim([1e-25, 5e-23])

    ax2.tick_params(
        left=False, labelleft=False,
        bottom=False, labelbottom=False
    )

    ax2.yaxis.set_minor_locator(ticker.NullLocator())

    mark_inset(ax1, ax2, loc1=1, loc2=3, fc="none", ec='0.5')

    plt.savefig("graph.png", dpi=300)

    if show:
        plt.show()


def calculate_differences(filenames: [str], dir: str):
    e_values = []

    for filename in filenames:
        e_values.append(compile_e_values(filename))

    file = open(f'data/{dir}/111SPOW.MAT')
    e_extp_line = file.readlines().pop()
    file.close()

    e_extp = Decimal(generate_values(e_extp_line)[1])

    result = [[] for i in range(2)]

    for index, values in enumerate(e_values):
        for value in values:
            result[index].append([
                value[0],
                (e_extp - value[1]) / e_extp
            ])

    return result


# Does the bulk of the work. Reads main file and compiles data
# Return value: [ [number, E-E_extp, Omega]1, []2, ..., []n ]
# For example, 111S251's entry would be: [ 251, E_111S251 - E_extp_111S251, Omega_111S251 ]
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
    lines.pop()

    for line in lines:
        e_values.append(generate_values(line))

    for value in e_values:
        e_values_converted.append([int(value[0]), Decimal(value[1])])

    return e_values_converted

# Extracts name and energy value from a given line.
# If search is true, also calls get_omega to get the relevant omega value
def generate_values(line: str):
    return [line.split()[0], line.split()[1]]

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