# Script to specify sets of model parameters and launch vertex model
# simulations of these sets as separate processes
import numpy as np

def calc_omega(taue, taul):
    """Calculate the predicted angular frequency of oscillation

    Prediction from linear stability of the 1d model presented in
    Boocock et al NatPhys (2021) DOI:10.1038/s41567-020-01037-7.

    taue: timescale for decay of ERK activation
    taul: timescale for decay of rest length / rest area

    Return: (double) angular frequency

    """
    p1 = 1/(taul**(1/2) * taue**(3/2))
    p2 = 1/(taul**(3/2) * taue**(1/2))
    p3 = 1/(taul*taue)
    return np.sqrt(p1 + p2 + p3)


def calc_period(taue, taul):
    """Calculate the predicted period of oscillation

    Prediction from linear stability of the 1d model presented in
    Boocock et al NatPhys (2021) DOI:10.1038/s41567-020-01037-7.

    taue: timescale for decay of ERK activation
    taul: timescale for decay of rest length / rest area

    Return: (double) period
    """
    return 2*np.pi / calc_omega(taue, taul)


if __name__ == "__main__":
    from run_batch import create_args_arr, launch_subprocess
    from multiprocessing import Pool
    from pathlib import Path

    # Path to executable within chaste_build (You may need to change
    # this line to specify the path to the executable on your own
    # system)
    executable1 = Path(*Path.cwd().parts[:-4],
                      "chaste_build/projects/erk_waves/test",
                      "TestERKWaveWithSelfPropulsionNoAlignment")

    # Specify a root directory (or stem) for simulation output. This
    # will appear within the directory specified by the environment
    # variable "CHASTE_TEST_OUTPUT" and contain output from the
    # simulations launched in this script.
    outdir_stem = "erk_waves_output/"

    # Specify the parameters used to label the output directory (a
    # separate file "params.txt" containing all parameters and values
    # will be generated withing the directory)
    outdir_labels = ["nx", "dt", "P0", "n_abcrit", "ab_ratio", "F0", "tau_p", "seed"]

    # Specify sets of model simulation parameters

    dt = 0.01    # timestep for vertex and ODE solver (I keep these the same)

    # Specify the parameters and values that are common across all
    # simulations
    const_arg_dict = {
        "-nx": 44,    # Number of cells in x
        "-ny": 44,    # Number of cells in y
        "-dt": dt,
        "-dt_ode": dt,

        # "-sampling_timestep_multiple": 5*int(0.1/dt),    # Save n*(ten times per T)

        # Initial mean values of variables
        "-init_erk": 0.0,
        "-init_A0": 1.0,    # Initial rest length in units of average cell area

        # Specify the noise magnitude (std) of initial perturbation in
        # vertex positions and ERK activity
        "-init_noise": 0.01,

        # Timescale associated with changes in preferred perimeter
        "-taul": 20.0,    # Non-dimensionalized by tauE=6min (timescale of ERK)
        # Persistence timescale of persistent random walk
        "-tau_p": 5.0,    # non-dimensionalized by tauE (30mins if tauE=6min)

        # Area and perimeter elasticity. Ratio of KA/KP taken from
        # Henkes et al NatCom (2020) (same as used in Saraswathibhatla
        # et al (2021)) with magnitudes (equiv substrate friction)
        # rescaled to give a wavelength of approx 20 cells. See table
        # in supplementary material of Boocock et al (2023)
        # doi:10.1101/2023.03.24.534111. The mechanical timescale
        # tau_r from the 1D model in Boocock et al (2021) goes like
        # ~1/(2*KA).
        "-KA": 4*0.3,
        "-KP": 4*0.09561,

        # Specify whether to check for internal intersections (0 for
        # True, 1 for False)
        "-check_for_internal_intersections": 0,

        # Self-propulsion force magnitude (or rather F0/\zeta, i.e. scaled by friction)
        "-F0": 0.15,
        # Ratio of mechanochemical coupling strengths alpha/beta
        "-ab_ratio": 18.5,
        # Magnitude of mechanochemical coupling strength alpha*beta as
        # a multiple of the prediction for the onset of stability from
        # linear stability analysis of the 1D model (see Boocock et al
        # 2021)
        "-n_abcrit": 1.8,
    }


    # Timescale of ERK activation used to non-dimensionalize
    tau_e = 1.0
    # Specify the length of the burn-in period (in terms of the
    # predicted period of oscillation)
    const_arg_dict["-end_time"] = np.round(100*calc_period(tau_e, const_arg_dict["-taul"]))
    # Specify the length of the data capture period (in terms of the
    # predicted period of oscillation)
    const_arg_dict["-bonus_time"] = np.round(10*calc_period(tau_e, const_arg_dict["-taul"]))
    # Specify a sampling interval for the data capture period
    const_arg_dict["-sampling_timestep_multiple"] = 5*int(0.1/dt)

    # Parameter lists within this dictionary are paired by index,
    # i.e. p1[i] with p2[i]. Lists should be of the same length.
    paired_arg_dict = {
    }


    # Dimensionless shape index P0/Sqrt(A0) controlling tissue
    # rigidity (see e.g. Bi et al NatPhys (2015)
    # https://doi.org/10.1038/nphys3471)
    p0 = np.array([
        3.5,
        3.8,
        ])

    # Run simulation for all combinations of parameter values within
    # the lists in this dictionary
    comb_arg_dict = {
        "-seed": [1],    # Random seeds
        "-P0": p0,
    }

    # Create arrays of arguments for different simulation
    # parameterizations to be passed to the executable and run as
    # separate processes
    args_arr = create_args_arr(executable,
                               const_arg_dict,
                               paired_arg_dict,
                               comb_arg_dict,
                               outdir_labels,
                               dir_name=outdir_stem)

    # Launch the simulations as separate subprocesses
    n_cpus = 2    # Number of simultaneous processes to run, i.e
                  # number of cores to use with each parameterization
                  # run on a separate core.
    with Pool(n_cpus) as p:
        p.map(launch_subprocess, args_arr)
