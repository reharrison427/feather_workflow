'''
Clean interferometer-only spectral line mosaic data

Initial shallow clean is designed for data with bright sidelobe emission
and strong negative bowls: sidelobethreshold and lownoisethreshold are 
lowered from their recommended ACA values to create an intial clean mask.
Checking initial mask is recommended to ensure that sidelobe emission and 
noise was not included in clean mask. 

Files produced by initial clean and ms files are copied to a directory 
called "iter1_niter5000". 

Deep clean raises sidelobethreshold and lownoisethreshold to recommended 
ACA values and cleans to noise-based threshold. 

CASA log parser detects warning messages related to diverging or negative
residuals, negative psf, empty residuals, or stopping criteria other than
hitting the iteration or noise limit. If warnings are found for initial 
shallow clean, the deep clean will not be run. 
'''

import glob
import re
import os
from pathlib import Path
import shutil
import numpy as np

# Helper functions
def get_latest_casa_log():
    """
    Return the most recently modified CASA log file in the
    current directory, excluding *_warnings.log files.
    """

    log_files = [
        f for f in glob.glob("casa*.log")
        if not f.endswith("_warnings.log")
    ]

    if not log_files:
        raise FileNotFoundError(
            "No casa*.log files found in current directory."
        )

    return max(log_files, key=os.path.getmtime)

def get_log_line_count(logfile):
    """
    Return the current number of lines in a CASA log file.
    Used to make sure parser only looks for warning messages after
    the last tclean started. 
    """

    with open(logfile, "r", errors="replace") as f:
        return sum(1 for _ in f)


def check_casa_log_since(logfile, start_line, stage):
    """
    Search only the portion of a CASA log written after
    start_line.

    If warnings are found, append them to the warning log under
    a label identifying the tclean stage.

    Line numbers refer to lines in the full CASA log.

    Parameters
    ----------
    logfile : str
        CASA log file to search.

    start_line : int
        Number of lines in the CASA log immediately before
        tclean was run.

    stage : str
        Label for this tclean stage, e.g. "INITIAL CLEAN"
        or "DEEP CLEAN".

    Returns
    -------
    found_warnings : bool
        True if any tclean warnings were found.
    """

    root, ext = os.path.splitext(logfile)
    warning_file = root + "_warnings.log"

    diverging_string = (
        "Exited multiscale minor cycle without reaching "
        "any stopping criterion"
    )

    negative_psf_string = (
        "As Peak of PSF is negative, you should setscales again "
        "with a smaller scale size"
    )

    global_stop_string = "Reached global stopping criterion"

    acceptable_stop_strings = (
        "iteration limit",
        "any combination of n-sigma and other valid exit criterion",
    )

    empty_residuals_string = ("Restoring with an empty model image")

    large_scale_divergence_string = ("Diverging due to large scale")

    channel_pattern = re.compile(r":C(\d+)\]")

    warnings = []

    # --------------------------------------------------------
    # Read CASA log
    # --------------------------------------------------------

    with open(logfile, "r", errors="replace") as f:

        for line_number, line in enumerate(f, start=1):

            # Ignore everything written before this tclean
            if line_number <= start_line:
                continue

            # ------------------------------------------------
            # Multiscale minor-cycle problem
            # ------------------------------------------------

            if diverging_string in line:

                match = channel_pattern.search(line)

                if match:
                    channel = int(match.group(1))

                    warnings.append(
                        f"Line {line_number}: "
                        f"WARNING: negative or diverging residuals "
                        f"in channel {channel}. "
                        f"Try decreasing cycleniter or increasing "
                        f"cyclefactor."
                    )

                else:
                    warnings.append(
                        f"Line {line_number}: "
                        f"WARNING: multiscale minor cycle exited "
                        f"without reaching a stopping criterion, "
                        f"but channel number could not be determined."
                    )

            # ------------------------------------------------
            # Negative PSF
            # ------------------------------------------------

            if negative_psf_string in line:

                warnings.append(
                    f"Line {line_number}: {line.rstrip()}"
                )

            # ------------------------------------------------
            # Unexpected global stopping criterion
            # ------------------------------------------------

            if global_stop_string in line:

                acceptable_stop = any(
                    text in line
                    for text in acceptable_stop_strings
                )

                if not acceptable_stop:

                    warnings.append(
                        f"Line {line_number}: "
                        f"WARNING: tclean stopped without reaching "
                        f"iteration or noise limit"
                    )

                    warnings.append(
                        f"Line {line_number}: {line.rstrip()}"
                    )

            # ------------------------------------------------
            # Empty residuals 
            # ------------------------------------------------

            if empty_residuals_string in line:

                warnings.append(
                    f"Line {line_number}: "
                    f"No emission in clean mask, try lowering automasking thresholds."
                )

                warnings.append(
                    f"Line {line_number}: {line.rstrip()}"
                )

            # ------------------------------------------------
            # Diverging due to large scale
            # ------------------------------------------------

            if large_scale_divergence_string in line:

                warnings.append(f"Line {line_number}: {line.rstrip()}")
                             

    # --------------------------------------------------------
    # Only create warning file if warnings were found
    # --------------------------------------------------------

    if warnings:

        file_exists = os.path.exists(warning_file)

        with open(warning_file, "a") as f:

            if not file_exists:
                f.write(
                    "Found tclean warning(s) in log file. "
                    "Check tclean parameters before proceeding.\n"
                )

            f.write(f"\n=== {stage} ===\n")

            for warning in warnings:
                f.write(warning + "\n")

        print(f"Saved warnings to: {warning_file}")

        return True

    return False

# Cleaning

'''
Initial shallow clean: low niter and loose automasking criteria
If noise limit is hit before iteration limit, lower niter.
If mask contains sidelobe emission, raise sidelobethreshold.
If mask contains noise, raise lownoisethreshold. 
'''
#######################################
# Check CASA version

import os
import glob
##################################################

# Image Parameters

# source parameters
# ------------------

field='3~40' # science field(s). For a mosaic, select all mosaic fields. DO NOT LEAVE BLANK ('') OR YOU WILL POTENTIALLY TRIGGER A BUG IN CLEAN THAT WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
# gridder='standard' # uncomment if single field 
gridder='mosaic' # uncomment if mosaic or if combining one 7m and one 12m pointing.
phasecenter=17


# image parameters.
# ----------------

cell='0.9arcsec' # cell size for imaging.
imsize = [300,230] # size of image in pixels.

# velocity parameters
# -------------------

outframe='lsrk' # velocity reference frame. 
veltype='radio' # velocity type. 

# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean. 

weighting = 'briggs'
robust=0.5


pbmask=0.0
noisethreshold=4.0
nsigma=3.0
#sidelobethreshold = 1.25
sidelobethreshold=0.9
lownoisethreshold = 1.5
minbeamfrac = 0.3
negativethreshold = 0.0
fastnoise = False
growiterations=75
dogrowprune=True

restfreq='230.538GHz' 

spw='0' # uncomment and replace with appropriate spw 

lineimagename =  'mosaic1_M17_CO_cube_I_int_cyclefactor1_scales05'

start='-17.0km/s' # start velocity. See science goals for appropriate value.
width='0.1835042268189km/s' # velocity width. See science goals.
nchan = 589
threshold='0.0Jy/beam'
niter=5000
CO_vis = glob.glob('*ms*.CO')

for linevis in CO_vis:
    clearcal(vis=linevis)
    delmod(vis=linevis)



tclean(vis=CO_vis,
       imagename=lineimagename,
       field=field,
       phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
       mosweight = True, # uncomment if mosaic
       specmode='cube',
       start=start,
       width=width,
       nchan=nchan,
       deconvolver='multiscale',
       scales=[0, 5], # scales in pixels.
       # Uncomment the below to image with nterms>1 when the fractional bandwidth is greater than 10%.
       #deconvolver='mtmfs',
       #nterms=2,
       imsize = imsize,
       cell= cell,
       weighting = weighting,
       robust = robust,
       niter = niter,
       threshold = threshold,
       interactive = False,
       gridder = 'mosaic',
       pbcor = True,
       pblimit=0.2,
       sidelobethreshold=sidelobethreshold,
       lownoisethreshold=lownoisethreshold,
       minbeamfrac=minbeamfrac,
       negativethreshold=negativethreshold,
       fastnoise=fastnoise,
       usemask="auto-multithresh",
       pbmask=pbmask,
       nsigma=nsigma,
       verbose=True,
       restoringbeam='common',
       parallel=True,
       outframe='lsrk', # velocity reference frame. 
       veltype='radio',
       restfreq=restfreq,
       restart=True)

'''
Copy files and ms's from inital clean to iter1_niter. 
'''
Path("iter1_niter5000").mkdir(exist_ok=True)
shutil.copytree("mosaic1_M17_CO_cube_I_int_cyclefactor1_scales05*", "iter1_niter5000", dirs_exist_ok=True)
shutil.copytree("*.ms*.CO", "iter1_niter5000", dirs_exist_ok=True)

# Check CASA log for warnings. If warnings are found, do not run deep clean.
casa_log = get_latest_casa_log()
line_count = get_log_line_count(casa_log)   
check_casa_log_since(casa_log, line_count, "INITIAL CLEAN")

if not check_casa_log_since(casa_log, line_count, "INITIAL CLEAN"):
    print("No warnings found in initial clean. Proceeding to deep clean.")
    '''
    Deep clean: clean down to noise-based threshold. 
    May need to be re-run if iteration limit is hit before
    noise-based threshold. 
    '''
    #######################################
    # Check CASA version

    import os
    import glob
    ##################################################

    # #############################################
    # Image Parameters


    # source parameters
    # ------------------

    field='3~40' # science field(s). For a mosaic, select all mosaic fields. DO NOT LEAVE BLANK ('') OR YOU WILL POTENTIALLY TRIGGER A BUG IN CLEAN THAT WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
    # gridder='standard' # uncomment if single field 
    gridder='mosaic' # uncomment if mosaic or if combining one 7m and one 12m pointing.
    phasecenter=17


    # image parameters.
    # ----------------

    cell='0.9arcsec' # cell size for imaging.
    imsize = [300,230] # size of image in pixels.

    # velocity parameters
    # -------------------

    outframe='lsrk' # velocity reference frame. 
    veltype='radio' # velocity type. 

    # imaging control
    # ----------------

    # The cleaning below is done interactively, so niter and threshold can
    # be controlled within clean. 

    weighting = 'briggs'
    robust=0.5

    pbmask=0.0
    noisethreshold=5.0
    nsigma=3.0
    sidelobethreshold = 1.25
    lownoisethreshold = 2.0
    minbeamfrac = 0.3
    negativethreshold = 0.0
    fastnoise = False
    growiterations=75
    dogrowprune=True
    restfreq='219.560358GHz'
    cycleniter=500
    spw='0' # uncomment and replace with appropriate spw 


    restfreq='230.538GHz' 

    spw='0' # uncomment and replace with appropriate spw 

    lineimagename =  'mosaic1_M17_CO_cube_I_int_cyclefactor1_scales05'

    start='-17.0km/s' # start velocity. See science goals for appropriate value.
    width='0.1835042268189km/s' # velocity width. See science goals.
    nchan = 589
    threshold='0.0Jy/beam'
    niter=100000
    CO_vis = glob.glob('*ms*.CO')

    tclean(vis=CO_vis,
        imagename=lineimagename,
        field=field,
        phasecenter=phasecenter, # uncomment if mosaic or imaging an ephemeris object
        mosweight = True, # uncomment if mosaic
        specmode='cube',
        start=start,
        width=width,
        nchan=nchan,
        deconvolver='multiscale',
        scales=[0, 5], # scales in pixels.
        # Uncomment the below to image with nterms>1 when the fractional bandwidth is greater than 10%.
        #deconvolver='mtmfs',
        #nterms=2,
        imsize = imsize,
        cell= cell,
        weighting = weighting,
        robust = robust,
        niter = niter,
        threshold = threshold,
        interactive = False,
        gridder = 'mosaic',
        pbcor = True,
        pblimit=0.2,
        sidelobethreshold=sidelobethreshold,
        lownoisethreshold=lownoisethreshold,
        minbeamfrac=minbeamfrac,
        negativethreshold=negativethreshold,
        fastnoise=fastnoise,
        usemask="auto-multithresh",
        pbmask=pbmask,
        nsigma=nsigma,
        verbose=True,
        restoringbeam='common',
        parallel=True,
        outframe='lsrk', # velocity reference frame. 
        veltype='radio',
        restfreq=restfreq,
        restart=True,
        cycleniter=cycleniter)

else:
    print("Warnings found in initial clean. Deep clean will not be run.")
    exit()



