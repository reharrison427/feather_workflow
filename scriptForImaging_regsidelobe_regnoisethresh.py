#######################################
# Check CASA version
'''
import re
import casadef

if casadef.casa_version < '4.4.0' :
    sys.exit("Please use CASA version greater than or equal to 4.4.0 with this script")
'''
import os
import glob
##################################################
# Create an Averaged Continuum MS


finalvis='calibrated_final.ms' # This is your output ms from the data
                               # preparation script.
'''
# Use plotms to identify line and continuum spectral windows.
plotms(vis=finalvis, xaxis='channel', yaxis='amplitude',
       ydatacolumn='data',
       avgtime='1e8', avgscan=True, avgchannel='1', 
       iteraxis='spw' )
'''

# Set spws to be used to form continuum
contspws = '0, 1, 2, 3, 4, 5'

# If you have complex line emission and no dedicated continuum
# windows, you will need to flag the line channels prior to averaging.
'''
flagmanager(vis=finalvis,mode='save',
            versionname='before_cont_flags')

initweights(vis=finalvis,wtmode='weight',dowtsp=True)

# Flag the "line channels"
flagchannels='0:230~600,2:10~25;54~60,3:45~53,4:450~580,5:500~565' 

flagdata(vis=finalvis,mode='manual',
          spw=flagchannels,flagbackup=False)

# check that flags are as expected, NOTE must check reload on plotms
# gui if its still open.
plotms(vis=finalvis,yaxis='amp',xaxis='channel',
       avgchannel='1',avgtime='1e8',avgscan=True,iteraxis='spw') 
'''
# Average the channels within spws
contvis='calibrated_final_cont.ms'
'''
rmtables(contvis)
os.system('rm -rf ' + contvis + '.flagversions')

split(vis=finalvis,
     spw=contspws,      
     outputvis=contvis,
      width=[1,1,8,8,1,1], # number of channels to average together. The final channel width should be less than 125MHz in Bands 3, 4, and 6 and 250MHz in Band 7.
     datacolumn='data')


# Check the weights. You will need to change antenna and field to
# appropriate values
plotms(vis=contvis, yaxis='wtsp',xaxis='freq',spw='',antenna='DA42',field='0')

# If you flagged any line channels, restore the previous flags
flagmanager(vis=finalvis,mode='restore',
            versionname='before_cont_flags')

# Inspect continuum for any problems
plotms(vis=contvis,xaxis='uvdist',yaxis='amp',coloraxis='spw')
'''
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
#sidelobethreshold=0.9
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
'''
for linevis in CO_vis:
    clearcal(vis=linevis)
    delmod(vis=linevis)
'''    
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




