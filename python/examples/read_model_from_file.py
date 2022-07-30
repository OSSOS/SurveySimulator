"""
Simulation to expose L7 model to CFEPS observation characterization
"""
import argparse
import sys

import ossssim

default_model_filename = '../../Models/L7model-3.0-0.0'
default_characterization_dir = '../../Surveys/CFEPS'


def run(model_filename: str = default_model_filename,
        characterization_dir: str = default_characterization_dir,
        seed: int = 113715191,
        track_condition: int = 0,
        output: str = 'detect.dat') -> None:
    """
    Run the survey simulator on these inputs.
    Args:
        model_filename: name of the SSurvey Simulator formatted model file, but needs to have a commented column header.
        characterization_dir: directory that contains the survey characterization pointing.list file
        seed: random number seed for the simulator to provide reproducible results.
        track_condition: how many source to track (-ve) or iteration to run (+ve) or do the entire file (0)
        output: name of file to store the detected sources to.
    Returns:
        None
    """
    model = ossssim.ModelFile(model_filename)
    ssim = ossssim.OSSSSim(characterization_dir)

    # create a detection output file and fill header with circumstance values that match the input model.
    detect_file = ossssim.DetectFile(output)
    detect_file.epoch = model.epoch
    detect_file.longitude_neptune = model.longitude_neptune
    detect_file.colors = model.colors
    detect_file.write_header(seed)

    # setup some variables to count the number of object tracked.
    n_iter = n_hits = n_track = 0

    # set track condition to stop when we have 1 tracked object
    # set value to -1 to run until 1 iterations has occurred
    # set value to 0 to loop to the end of the file.
    for row in model:
        n_iter += 1
        result = ssim.simulate(row, seed=seed, epoch=model.epoch)
        if result['flag'] > 0:
            n_hits += 1
            detect_file.write_row(result)
        # stop the loop when maximum/desired detections have occurred if ntrack==0 then go to end of model file.
        if (0 < track_condition <= n_track) | (0 < -track_condition <= n_iter):
            break
    detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)


def main():
    """
    The main console script.  Read SSim input parameters from command-line arguments or from a Drive.in file.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="This script can be called with arguments or can read Driver.in from standard input.",
        epilog='''\n\n---- OR READ FROM STDIN -----
Example of Driver.in format:

123456789                    ! Random number generator seed
0                            ! Number of detections to track or draws to take (0 == do the entire file)
Surveys/CFEPS                ! Survey description directory name
Models/L7model-3.0-9.0       ! Name of model parameter file
./L7Detected.dat             ! Output file name for all simulated detections
./L7Tracked.dat              ! Output file name for characterized, tracked simulated detections
               
               
USAGE when reading from standard input: read_model_from_file < Driver.in''')
    parser.add_argument('model',
                        nargs='?', default=default_model_filename,
                        type=str,
                        help="File with SSim input model, should have a commented header line")
    parser.add_argument('characterization', nargs='?',
                        default=default_characterization_dir,
                        type=str,
                        help="Directory with the characterization files")
    parser.add_argument('--seed', type=int, default=113715191, help='Random seed to pass to Survey Simulator')
    parser.add_argument('--track_condition', type=int, default=0,
                        help="If +ve then run until that number of objects tracked, "
                             "if -ve then run for that many iterations, if 0 then run to end of model")
    parser.add_argument('--output', type=str, help="Name of file to store detected/tracked sources to.")

    if len(sys.argv) > 1:
        args = parser.parse_args()
        model_filename = args.model
        characterization_dir = args.characterization
        seed = args.seed
        track_condition = args.track_condition
        output_filename = args.output
    else:
        seed = int(sys.stdin.readline().split()[0])
        # set track condition to +ve value to stop when we have +ve tracked object
        # set value to -ve value to run until that number of iterations has occurred
        # set value to 0 to loop to the end of the file.
        track_condition = int(sys.stdin.readline().split()[0])
        characterization_dir = sys.stdin.readline().split()[0]
        model_filename = sys.stdin.readline().split()[0]
        output_filename = sys.stdin.readline().split()[0]

    run(model_filename=model_filename,
        characterization_dir=characterization_dir,
        seed=seed,
        track_condition=track_condition,
        output=output_filename)


if __name__ == '__main__':
    main()
