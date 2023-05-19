# Functions to create argument arrays of parameter values and launch
# simulations as separate subprocesses
import itertools as it
import subprocess


def create_args_arr(executable, const_arg_dict, paired_arg_dict,
                    comb_arg_dict, outdir_labels=None, dir_name=""):
    """Create an array of different simulation configurations

    output: list - Each row is a set of arguments (execuateble and set
    of options) to provide to subprocess.run() and launch a different
    simulation.
    """
    # Add constant parameters
    const_args = [(k, str(v)) for k, v in const_arg_dict.items()]

    # Add pairs of variable parameters
    args_arr = []
    try:
        n = len(next(iter(paired_arg_dict.values())))
        for i in range(n):
            args_arr.append(const_args + [(k, str(v[i])) for k, v in paired_arg_dict.items()])
    except StopIteration:
        pass

    if not args_arr:
        args_arr = [const_args]

    # Create combinations of variable parameters
    for k, vals in comb_arg_dict.items():
        args_arr = [row + [(k, str(v))] for v in vals for row in args_arr]

    # Add an option for the output directory name
    if outdir_labels is not None:
        # Name directories as "label1_value1_label2_value2..."
        for i, row in enumerate(args_arr):
            dname = dir_name
            for label in outdir_labels:
                d = dict(row)
                if label == "ab_ratio":
                    # Deal with long trailing numbers
                    label = label + "_" + "{0:1.4f}".format(float(dict(row)["-" + label]))
                elif label == "P0":
                    # Deal with long trailing numbers
                    label = label + "_" + "{0:1.4f}".format(float(dict(row)["-" + label]))
                else:
                    label = label + "_" + str(dict(row)["-" + label])
                if len(dname) == 0:
                    dname += label
                else:
                    dname += "_" + label
            args_arr[i].append( ("-outdir", dname) )
    else:
        # Label output directories with unique integers
        n_didgets = len(str(len(args_arr)))
        fmt_str = "simulation_output_{0:0" + str(n_didgets) + "}"
        for i, _ in enumerate(args_arr):
            args_arr[i].append( ("-outdir", fmt_str.format(i)) )

    # Add the executable to the start and combine each row into one array
    args_arr = [[executable] + list(it.chain.from_iterable(row)) for row in args_arr]
    return args_arr


def launch_subprocess(args):
    try:
        p = subprocess.run(args=args, check=True)
    except Exception:
        print("Process not complete ", args)
