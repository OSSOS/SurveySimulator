""""
This example uses the ssim_util functions to open a model file, passes through detect and then writes the result to a detect file.

Run this as a script via command:

python read_model_from_file.py < ReadModelFromFile.in
"""
import numpy
import sys
import osssim

def run(model_filename,
        detect_filename,
        track_filename,
        characterization_directory,
        seed,
        track_condition):
    """
    Loop over the model_file until we hve track detections.
    """
    model = osssim.SSimModelFile(model_filename, randomize=False)
    ssim = osssim.OSSSSim(characterization_directory)

    # create a detection output file and fill header with circumstance values that match the input model.
    detect_file = osssim.DetectFile(detect_filename)
    detect_file.epoch = model.epoch.jd
    detect_file.lambda_neptune = model.longitude_neptune
    detect_file.colors = model.colors
    detect_file.write_header(seed)

    track_file = osssim.DetectFile(track_filename)
    track_file.epoch = model.epoch.jd
    track_file.lambda_neptune = model.longitude_neptune
    track_file.colors = model.colors
    track_file.write_header(seed)


    # setup some variables to count the number of object tracked.
    n_iter = n_hits = n_track = 0

    # set track condition to stop when we have 1 tracked object
    # set value to -1 to run until 1 iterations has occurred
    # set value to 0 to loop to the end of the file.
    for row in model:
        n_iter += 1
        result = ssim.simulate(row, seed=seed, epoch=model.epoch.jd)
        if result['flag'] > 0:
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


if __name__ == '__main__':
    seed = int(sys.stdin.readline().split()[0])
    track_condition = int(sys.stdin.readline().split()[0])
    characterization_dir = sys.stdin.readline().split()[0]
    model_filename = sys.stdin.readline().split()[0]
    detect_filename = sys.stdin.readline().split()[0]
    track_filename = sys.stdin.readline().split()[0]
    run(model_filename, detect_filename, track_filename,
        characterization_dir, seed, track_condition)
