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
