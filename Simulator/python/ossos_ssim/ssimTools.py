"""ssimTools is a module called by Driver.py. ssimTools contains the miscellaneous functions required to run Driver,
primarily the file input/output functions. """
import numpy as np


def input_variables(f_name):
    """
    Opens and reads in the required input files before returning them to Driver.

    Parameters
    ----------
    f_name - the name of the input file to be read in.

    Return Values
    -------------
    obj_file - The name of the arbitrarily named module of the form of the GiMeObj module. (no file extension)
    seed - the number used to seed the random number generator in GiMeObj.
    n_detect_max - Number of objects to detect before Driver will terminate.
    n_track_max - Number of objects to track before Driver will terminate.
    survey_dir -  Path to directory containing the survey characterizations.
    """
    with open(f_name, 'r') as f:
        tmp = f.read().split()
        obj_file = tmp[0]
        seed = int(tmp[1])
        n_detect_max = int(tmp[2])
        n_track_max = int(tmp[3])
        survey_dir = tmp[4]

    # This section can be used to draw the input variables from several different files. Useful on a cluster.
    ####################################################################################
    # with open('input.file', 'r') as f:
    #     tmp = f.read().split()
    #     obj_file = tmp[0]  # Name of GiMeObj module file (without .py file extension)
    #
    # with open('seeds.txt', 'r') as f:
    #     seed = int(f.read())  # Seed for random number generator
    #
    # with open('number_to_detect.txt', 'r') as f:
    #     n_detect_max = int(f.readline())  # Number of objects to detect
    #
    # with open('number_to_track.txt', 'r') as f:
    #     n_track_max = int(f.readline())  # Number of objects to track
    #
    # with open('surveydir.txt', 'r') as f:
    #     survey_dir = f.read().split()[0]  # Path to directory containing the characterization files
    ####################################################################################

    return obj_file, seed, n_detect_max, n_track_max, survey_dir


# Functions to create output files and draw their headers.


def drawn_header(f_name):
    """
    Open drawn file and write header

    The file format of the header set by SurveySimulator.
    """
    # TODO make this header for the columns that will be written.
    with open(f_name, 'w') as f:
        f.write(f'# {"a":>7} {"e":>6} '
                f'{"inc":>8} {"node":>8} '
                f'{"peri":>8} {"Manom":>8} '
                f'{"H":>6} {"resamp":>8} '
                f'{"Comments":>10s}\n'
                )


def det_header(f_name, seed):
    """
    Open detection file and write header
    """
    with open(f_name, 'w') as f:
        f.write(
            f'# Seed: {seed:10d}\n#\n'
            f'# flag: >0: detected; >2: characterized; 0 mod(2): tracked\n#\n'
            f'#{"a [AU]":>7s} {"e":>6s} {"i [°]":>8s} {"Ω [°]":>8s} {"ω [°]":>8s} {"M [°]":>8s} '
            f'{"ResAmp [°]":>10s} {"q [°]":>8s} {"Dist_* [AU]":>12s} {"M(t) [°]":>8s} {"MagR":>8s} '
            f'{"H MagR":>6s} {"Color":>5s} {"flag":>5s} {"Dist_E [AU]":>12s} {"Mag_Intr":>8s} '
            f'{"H_Intr":>6s} {"eff":>4s} {"RA [H]":>8s} {"Dec [°]":>8s} \u202F{"𝛿RA [″/hr]":>11s} '
            f'\u200A{"𝛿Dec [″/hr]":>11s} {"Survey":>10s} {"Comments":>10s}\n#\n'
        )
        # Special Unicode spaces (U+202F Narrow No-Break Space, and U+2009 Thin Space) are needed to ensure proper
        # spacing in the text file. This is because the double prime symbol used to denote arcseconds is thinner than
        # standard fixed width unicode characters.


def track_header(f_name):
    """
    Open tracked detection file and write header
    """
    with open(f_name, 'w') as f:
        f.write(
            f'#{"a [AU]":>7s} {"e":>6s} {"i [°]":>8s} {"Ω [°]":>8s} {"ω [°]":>8s} {"M [°]":>8s} '
            f'{"ResAmp [°]":>10s} {"q [°]":>8s} {"Dist_* [AU]":>12s} {"M(t) [°]":>8s} {"MagR":>8s} '
            f'{"H MagR":>6s} {"Color":>5s} {"Comments":>10s}\n#\n'
        )


# Functions to append and write data to the already created files with headers.


def drawn_write(f_name, a, e, inc, node, peri, M, h, resamp, comments):
    """
    Append data to the drawn file.
    """
    with open(f_name, 'a') as f:
        f.write(
            f"{a:8.3f} {e:6.3f} {np.rad2deg(inc):8.3f} {np.rad2deg(node):8.1f} "
            f"{np.rad2deg(peri):8.1f} {np.rad2deg(M):8.1f} {h:6.2f} {np.rad2deg(resamp):8.1f} "
            f"{comments:10s}"
            f"\n"
        )


def det_write(f_name, a, e, inc, node, peri, m, resamp, r, mt, m_rand, h_rand, color, ic, flag,
              delta, m_int, h, eff, ra, dec, d_ra, d_dec, surna, comments):
    """
    Append data to the detections file.
    """
    with open(f_name, 'a') as f:
        f.write(
            f'{a:8.3f} {e:6.3f} {np.rad2deg(inc):8.3f} {np.rad2deg(node):8.1f} {np.rad2deg(peri):8.1f} '
            f'{np.rad2deg(m):8.1f} {np.rad2deg(resamp):10.5f} {a * (1 - e):8.3f} {r:12.5f} {np.rad2deg(mt):8.3f} '
            f'{m_rand:8.3f} {h_rand:6.2f} {color[ic - 1]:5.2f} {flag:>5d} {delta:12.5f} {m_int:8.3f} '
            f'{h:6.2f} {eff:4.2f} {np.rad2deg(ra) / 15:8.5f} {np.rad2deg(dec):8.4f} '
            f'{np.rad2deg(d_ra) * 3600 / 24:11.6f} {np.rad2deg(d_dec) * 3600 / 24:11.6f} '
            f'{surna.decode("utf-8"):>10s} {comments:>10s}\n'
        )


def track_write(f_name, a, e, inc, node, peri, m, resamp, r, mt, m_rand, h_rand, color, ic, comments):
    """
    Append data to the tracked file.
    """
    with open(f_name, 'a') as f:
        f.write(
            f'{a:8.3f} {e:6.3f} {np.rad2deg(inc):8.3f} {np.rad2deg(node):8.1f} {np.rad2deg(peri):8.1f} '
            f'{np.rad2deg(m):8.1f} {np.rad2deg(resamp):10.5f} {a * (1 - e):8.3f} {r:12.5f} {np.rad2deg(mt):8.3f} '
            f'{m_rand:8.3f} {h_rand:6.2f} {color[ic - 1]:5.2f} {comments:>10s}\n'
        )


# Function to append a suffix when Driver terminates.


def det_suffix(f_name, n_drawn, n_det, n_track):
    """
    Append suffix data to the end of the detections file. It lists the total numbers of objects read-in, detected,
    and tracked, respectively.
    """
    with open(f_name, 'a') as f:
        f.write(
            '#\n'
            f'# Total number of objects:    {n_drawn:11d}\n'
            f'# Number of detected objects: {n_det:11d}\n'
            f'# Number of tracked objects:  {n_track:11d}\n'
        )
