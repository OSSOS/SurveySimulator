""""
This example runs the simulator from a driver file, as is done in the fortran usage.

Run this as a script via command:

python read_model_from_file.py < Driver.in
"""
import numpy
import sys
import ossssim


def run(model_filename,
        detect_filename,
        track_filename,
        characterization_directory,
        seed,
        track_condition):
    """
    Loop over the model_file until we have achieved the track_condition.
    """
    model = ossssim.ModelFile(model_filename)
    model.colnames.append('j')
    model.colnames.append('k')
    ssim = ossssim.OSSSSim(characterization_directory)

    # create a detection output file and fill header with circumstance values that match the input model.
    detect_file = ossssim.DetectFile(detect_filename)
    detect_file.epoch = model.epoch
    detect_file.lambda_neptune = model.longitude_neptune
    detect_file.colors = model.colors
    detect_file.write_header(seed)

    track_file = ossssim.DetectFile(track_filename)
    track_file.epoch = model.epoch
    track_file.lambda_neptune = model.longitude_neptune
    track_file.colors = model.colors
    track_file.write_header(seed)

    # setup some variables to count the number of object tracked.
    n_iter = n_hits = n_track = 0

    for row in model:
        n_iter += 1
        result = ssim.simulate(row, seed=seed, epoch=model.epoch)
        #if result['flag'] > 0:
        n_hits += 1
        detect_file.write_row(result)

        if (result['flag'] > 2) and (numpy.mod(result['flag'], 2) == 0):
            n_track += 1
            track_file.write_row(result)
        # stop the loop when maximum/desired detections have occurred if ntrack==0 then go to end of model file.
        if (0 < track_condition <= n_track) | (0 < -track_condition <= n_iter):
            break
    detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)
    track_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)


def parse_stdin():
    """
    Parse the stdin, expecting a Driver formatted file.

    Example of driver.in format:

    123456789		                   ! Random number generator seed
    0			                       ! Number of tracked detections
    Surveys/CFEPS                ! Survey description directory name
    Models/L7model-3.0-9.0       ! Name of model parameter file
    ./L7Detected.dat                  ! Output file name for all simulated detections
    ./L7Tracked.dat                   ! Output file name for characterized, tracked simulated detections

    """
    seed = int(sys.stdin.readline().split()[0])
    # set track condition to +ve value to stop when we have +ve tracked object
    # set value to -ve value to run until that number of iterations has occurred
    # set value to 0 to loop to the end of the file.
    track_condition = int(sys.stdin.readline().split()[0])
    characterization_dir = sys.stdin.readline().split()[0]
    model_filename = sys.stdin.readline().split()[0]
    detect_filename = sys.stdin.readline().split()[0]
    track_filename = sys.stdin.readline().split()[0]
    run(model_filename, detect_filename, track_filename,
        characterization_dir, seed, track_condition)


if __name__ == '__main__':
    parse_stdin()
