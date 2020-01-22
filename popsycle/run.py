#! /usr/bin/env python
"""
run.py
Run the entire PopSyCLE pipeline starting
"""

def write_galaxia_params(output_dir, output_basename='ZTF1',
                         seed=0,
                         longitude=45.19260648,
                         latitude=4.93717557,
                         area=10.00):
    params = [
        "outputFile %s" % output_basename,
        "outputDir %s" % output_dir,
        "photoSys UBV",
        "magcolorNames V,B-V",
        "appMagLimits[0] -1000",
        "appMagLimits[1] 1000",
        "absMagLimits[0] -1000",
        "absMagLimits[1] 1000",
        "colorLimits[0] -1000",
        "colorLimits[1] 1000",
        "geometryOption 1",
        "longitude %f" % longitude,
        "latitude %f" % latitude,
        "surveyArea %.2f" % area,
        "fSample 1",
        "popID -1",
        "warpFlareOn 1",
        "seed %i" % seed,
        "r_max 30",
        "starType 0",
        "photoError 0"
    ]
    with open(output_dir + '/galaxia_params.%i.txt' % seed, 'w') as f:
        for param in params:
            f.write(param + '\n')
