#! /usr/bin/env python3
import logging

import argparse
import sys
sys.path.append("../fortran")
import ssim_util
import SurveySubsF95
import numpy as np


def conv(st):
    return "".join(chr(x) for x in st)


SS = SurveySubsF95.surveysub
datadec = SurveySubsF95.datadec


# Open tracked detection file (no header)

def simulate(configs, gb=-0.12, ph=0, period=0.0, amp=0.0):
        """
        configs: a configuation object returned from ssim_util.SimConfigs
        gb: Bowell surge effect value
        ph: phase of rotation at epoch
        period: period of rotation
        amp: light-curve amplitude
        """

        detect_file = ssim_util.DetectFile(configs.detect_filename)
        detect_file.write_header(configs.seed)
        detect_file.write_header(configs.seed)
        # if ntrack is greater than 0 we are stopping at ntrack detections
        # if ntrack is less than 0 we stop of that many iterations.
        # for those case we randomize our model input
        model_table = ssim_util.SSimModelReader(configs.model_filename, randomize=(configs.ntrack != 0))
        n_iter = n_hits = n_track = 0

        for row in model_table:
            o_m = datadec.t_orb_m()
            for colname in row:
                if hasattr(o_m, colname):
                    setattr(o_m, colname, row[colname])
                    if colname in ['node', 'peri', 'm', 'inc']:
                        setattr(o_m, colname, np.radians(getattr(o_m, colname)))

            # attempt to detect this object.

            row['flag'], row['RA'], row['DEC'], row['d_ra'], row['d_dec'], row['r'], row['delta'],\
            row['m_int'], row['m_rand'], row['eff'], isur, row['M'], jdayp, ic, row['Survey'], row['h_rand'], ierr = \
                SS.detos1(o_m, model_table.keywords['epoch'],
                          row['h'], model_table.keywords['colors'],
                          gb, ph, period, amp,
                          configs.survey_dir, configs.seed)
            if ierr != 0:
                raise IOError("call to detos1 failed: {ierr}")
            # ic gives the filter that the object was 'detected' in, this allows us to determine the color
            if row['flag'] > 0:
                row['color'] = model_table.keywords['colors'][ic - 1]
                row['q'] = row['a'] * (1 - row['e'])
                row['M'] = np.degrees(row['M'])
                row['RA'] = np.degrees(row['RA'])/15.0
                row['DEC'] = np.degrees(row['DEC'])
                row['Survey'] = row['Survey'].decode('utf-8')
                # m_int and h are in "x" band (filter of object creation)
                # m_rand and h_rand are in discovery filter
                n_hits += 1
                detect_file.write_row(row)
                if (row['flag'] > 2) and (np.mod(row['flag'], 2) == 0):
                    n_track += 1

            # stop the loop when maximum/desired detections have occurred if configs.ntrack==0 then go to end of model file.
            if (0 < configs.ntrack <= n_track) | (0 < -configs.ntrack <= n_iter):
                break
        detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    logging.basicConfig(level=logging.DEBUG)
    parser.add_argument('config', help="Input configuration file.")
    args = parser.parse_args()
    configs = ssim_util.SSimConf(args.config)
    simulate(configs)

