#!/usr/bin/env python

"""
Observational SNR and exposure time calculator
Jonathan Aguilar
Nov. 15, 2013


license
----------
Copyright 2014 Jonathan A. Aguilar (JHU)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


usage
----------
- Run from ipython with %run ccdnoise.py or from command line with ./ccdnoise.py
- Click on the sliders to set the parameter values. 


to-do
----------
- Convert to a proper GUI backend to improve performance
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider#, Button, RadioButtons
import warnings

# some text formatting
mpl.rcParams["font.size"]=16
mpl.rcParams["text.color"]='white'
mpl.rcParams["axes.labelcolor"]='white'
mpl.rcParams["xtick.color"]='white'
mpl.rcParams["ytick.color"]='white'

def mag2flux(mag): 
    """
    For quickly converting visible-band magnitudes into photons/(um m^2 sec)
    """
    return 9.6e10*(10**(-mag/2.5))
def flux2mag(flux):
    """
    For quickly converting photons/(um m^2 sec) into V-band magnitudes
    """
    return -2.5*np.log10(flux/9.6e10)

def annotateAxis(xvals, yvals, labels, ax):
    """
    Semi-intelligently add labels to the plots without too much overcrowding
    """
    for label, x, y, i in zip( labels, xvals, yvals, range(targets.size) ):
        # check if label index is even or odd
        aboveBelow = i%2 and True or False # true if odd
        if aboveBelow:
            aboveBelow= -1
        else:
            aboveBelow = 1
        ax.annotate(
            label,
            color="black",
            xy=(x, y), xytext=(40, aboveBelow*10),
            textcoords='offset points', ha='right', va='bottom',
            #bbox=dict(boxstyle = 'round,pad=0.5', fc = 'yellow', fill=False, alpha = 0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))   


#plt.ion()

# PLOTTING
gs = gridspec.GridSpec(21,16,left=0.07,right=0.97)
fig = plt.figure()
fig.suptitle("SNR Calculator for CCDs")
fig.set_facecolor("black")

"""
# background image - disabled for now; was a cool image of a galaxy
bgndimg = mpl.image.imread("NGC1097HaLRGBpugh.jpg")
fig.figimage(bgndimg)
fullfigax= fig.add_axes([0,0,1,1])
fullfigax.imshow(bgndimg,)
"""

col1header = fig.text(0.2,0.92,r"Sources and noise",size="large")
col2header = fig.text(0.65,0.92,r"Telescope parameters",size="large")
# flux column width
col1width=7
col2width=4


# common x-axis

#xlims = (np.log10(1e-4), np.log10(1e15))
xlims = (np.log10(1e-2), np.log10(1e10))
middle = np.mean(xlims)

dt = np.dtype([("flux",'f8'),("label","S20")])

# Astronomical targets [Photons/(um m^2 sec)]
targets = np.array([
        (4.0e12, "Jupiter"),
        (0.1, "HST deep field"),
        (7.5, "Fomalhaut b"),
        #(4.0e9, "Andromeda"),
        (6.6e5, "QSO 3C 273"),
        (mag2flux(5.8),"GRB 080319B"),
        (mag2flux(24), "SN 1997ck"),
        (mag2flux(0), "Vega")
        ], dtype=dt)
targets.sort(order="flux")

targetsax = plt.subplot(gs[2:4,:col1width])#plt.sca(axes[numaxis_targets])
plt.xlim(xlims[0],xlims[1])
plt.ylabel("Targets",rotation=45)
plt.xlabel(r"$log_{10}[N_{\gamma}/ (m^2 \mu m\ sec)]$", labelpad=-1)
plt.scatter(np.log10(targets["flux"]),
                  0.5*np.ones(targets["flux"].size))
annotateAxis(np.log10(targets["flux"]),
             0.5*np.ones(targets["flux"].size),
             targets["label"], 
             plt.gca())
targetsax.yaxis.set_ticks([])
targetsax_mag = targetsax.twiny()
targetsax_mag.yaxis.set_ticks([])
targetsax_mag.set_xticklabels(['%1.1f'%flux2mag(10**(t)) for t in targetsax.get_xticks()])
targetsax_mag.set_xlim(xlims)
targetsax_mag.set_xlabel("mag", labelpad=1)
mytarget=np.array([(1e4,"Target")],dt)
mytargplt, = targetsax.plot(np.log10(mytarget["flux"][0]), 0.5, 
                            c='red',marker='x')
      
stargets=Slider(plt.subplot(gs[0,:col1width]),
                "",
                 xlims[0],xlims[1],
                 valinit=np.log10(mytarget["flux"][0]))


# Sky [Photons/(um m^2 sec)]
sky = np.array([
        (6.1e3,"Full moon"), 
#        (5.0e2, "Half moon"), 
        (2.9e2, "New moon"),
        ], dt)
sky.sort(order="flux")

skyax = plt.subplot(gs[7:9,:col1width])
plt.xlim(xlims[0],xlims[1])
plt.ylabel("Sky",rotation=45)
plt.xlabel(r"$log_{10}[N_{\gamma}/ (m^2 \mu m\ sec)]$", labelpad=-1)
plt.scatter(np.log10(sky["flux"]),
                  0.5*np.ones(sky["flux"].size))
annotateAxis(np.log10(sky["flux"]),
             0.5*np.ones(sky["flux"].size),
             sky["label"], 
             plt.gca())
skyax.yaxis.set_ticks([])
skyax_mag = skyax.twiny()
skyax_mag.yaxis.set_ticks([])
skyax_mag.set_xticklabels(['%1.1f'%flux2mag(10**(t)) for t in skyax.get_xticks()])
skyax_mag.set_xlim(xlims)
skyax_mag.set_xlabel("mag", labelpad=1)

mysky=np.array([(3e2,"Sky")],dt)
myskyplt, = skyax.plot(np.log10(mysky["flux"][0]), 0.5, 
                       c='red',marker='x')
   
ssky=Slider(plt.subplot(gs[5,:col1width]),
            "",
            xlims[0],xlims[1],
            valinit=np.log10(mysky["flux"][0]),
            facecolor="red")


# Dark current, [e-/(pixel sec)]
darkcurrent = np.array([
        (6.2e-3, "HST WFC"),
        (8e-1, "Keck NIRSPEC"),
        (1e-1, "Commercial CCD")
#        (1e-3,"40 C"), 
#        (1e-1,"80 C"), 
#        (1e2, "100 C")
        ], dt)
darkcurrent.sort(order="flux")

darkcurrentax = plt.subplot(gs[11:13,:col1width])
plt.xlim(-3,10)
plt.ylabel("Dark current",rotation=45)
plt.xlabel(r"$log_{10}[N_{e^-}/ (pixel \ sec)]$", labelpad=-1)
plt.scatter(np.log10(darkcurrent["flux"]),
            0.5*np.ones(darkcurrent["flux"].size))
annotateAxis(np.log10(darkcurrent["flux"]),
             0.5*np.ones(darkcurrent["flux"].size),
             darkcurrent["label"], 
             darkcurrentax)
darkcurrentax.yaxis.set_ticks([])
mydarkcurrent=np.array([(1e-2,"Dark current")],dt)
mydarkcurrentplt, = darkcurrentax.plot(np.log10(mydarkcurrent["flux"][0]), 0.5, 
                                       c='red',marker='x')
sdarkcurrent=Slider(plt.subplot(gs[10,:col1width]),
              "",
              darkcurrentax.get_xlim()[0],darkcurrentax.get_xlim()[1],
              valinit=np.log10(mydarkcurrent["flux"][0]),
              facecolor="red")


# Read noise, [e-/pix/read]
readnoise = np.array([
        (2., "Chandra ACIS"),
        (3., "HST WFC3"),
        (4., "P1640"),
        (10., "Subaru Suprime-Cam")
        ], dt)

readnoiseax = plt.subplot(gs[15:17,:col1width])
plt.xlim(0,15)
plt.ylabel("Read noise",rotation=45)
plt.xlabel(r"$N_{e^-}/pixel$", labelpad=-1)
plt.scatter(readnoise["flux"],
            0.5*np.ones(readnoise["flux"].size))
annotateAxis(readnoise["flux"],
             0.5*np.ones(readnoise["flux"].size),
             readnoise["label"],
             readnoiseax)
readnoiseax.yaxis.set_ticks([])
myreadnoise=np.array([(2,"Read noise")],dt)
myreadnoiseplt, = readnoiseax.plot(myreadnoise["flux"], 0.5,
                                   c="red",marker="x")
sreadnoise=Slider(plt.subplot(gs[14,:col1width]),
                  "",
                  readnoiseax.get_xlim()[0],readnoiseax.get_xlim()[1],
                  valinit=2,
                  facecolor="red")


# Cosmics, [N/(cm^2 sec)]
cosmics = np.array([
        (1.e-2,"Sea lvl"),
        (2e-2, "4 km"),
        (1,"Space"),
        (1.3e-4, "LHC")
        ], dt)
cosmics.sort(order="flux")

cosmicsax = plt.subplot(gs[19:21,:col1width])
plt.xlim(-5,1)
plt.ylabel("Cosmics",rotation=45)
plt.xlabel(r"$log_{10}[N_{cosmics}/ (cm^2\ sec)]$", labelpad=-1)
plt.scatter(np.log10(cosmics["flux"]),
            0.5*np.ones(cosmics["flux"].size))
annotateAxis(np.log10(cosmics["flux"]),
             0.5*np.ones(cosmics["flux"].size),
             cosmics["label"], 
             cosmicsax)
cosmicsax.yaxis.set_ticks([])
mycosmics=np.array([(1e-2,"Cosmics")],dt)
mycosmicsplt, = cosmicsax.plot(np.log10(mycosmics["flux"][0]), 0.5, 
                                       c='red',marker='x')
scosmics=Slider(plt.subplot(gs[18,:col1width]),
                "",
                cosmicsax.get_xlim()[0],cosmicsax.get_xlim()[1],
                valinit=np.log10(mycosmics["flux"][0]),
                facecolor="red")


# Telescope parameters
saperture=Slider(plt.subplot(gs[0,col1width+4:col1width+4+col2width]),
                 r"Aperture diameter [$m$]",
                 0,10,
                 valinit=1,
                 facecolor="green")
sbandwidth=Slider(plt.subplot(gs[1,col1width+4:col1width+4+col2width]),
                  r"Bandwidth [$\mu m$]",
                  0,10,
                  valinit=1,
                  facecolor="green")
sccdsize=Slider(plt.subplot(gs[2,col1width+4:col1width+4+col2width]),
                r"CCD dimension [pixels]",
                500,5000,
                valinit=2048,
                valfmt="%i",
                facecolor="green")
sqe=Slider(plt.subplot(gs[3,col1width+4:col1width+4+col2width]),
           r"Quantum efficiency",
           0,1,
           valinit=0.85,
           facecolor="green")
spixelsize=Slider(plt.subplot(gs[4,col1width+4:col1width+4+col2width]),
                  r"Pixel size [$\mu m$]",
                  1,100,
                  valinit=18,
                  facecolor="green")
spixelsperpsf=Slider(plt.subplot(gs[5,col1width+4:col1width+4+col2width]),
                     r"Pixels under PSF [$\#$]",
                     1,40,
                     valinit=10,
                     valfmt='%1.1f',
                     facecolor="green")
spsfsontarget=Slider(plt.subplot(gs[6,col1width+4:col1width+4+col2width]),
                     r"Target size [$PSFs$]",
                     1,200,
                     valinit=10,
                     facecolor="green")

# signal-to-noise
def calcSNR(sig=1,bgnd=0,readnoise=0,Idc=0,npix=1,time=1):
    """
    This is the famous CCD equation.
    The number of pixels is already taken into account 
    if you pull values from the bar chart, so leave it at 1
    """
    snr = []
    time = np.asarray(time)
    warnings.filterwarnings("error")
    try:        
        snr = [sig*t/np.sqrt( (sig+bgnd*npix+Idc*npix)*t+(readnoise**2)*npix ) for t in time]
    except RuntimeWarning:
        print "Error: cannot handle time = 0. Try again."
    return np.array(snr)

def calcCosmicHits(time=1):
    """
    The probability of losing a pixel to cosmics
    """
    time = np.asarray(time)
    ccdflux = (spixelsize.val*sccdsize.val*1e-4)**2 * 10**scosmics.val # flux through ccd
    npix = spsfsontarget.val*spixelsperpsf.val
    probPerPix = (1-1./(sccdsize.val**2))**(ccdflux*time)
    probPerTargetPix = probPerPix**npix
    return (probPerPix,probPerTargetPix)
    
# Signal size bar plot - the contribution of electrons from each source
signalsvals = np.array([
        ((10**stargets.val)*(np.pi*saperture.val**2)*sbandwidth.val*sqe.val,"Target"),
        ((10**ssky.val)*(np.pi*saperture.val**2)*sbandwidth.val*sqe.val,"Sky"),
        ((10**scosmics.val) * (spixelsize.val)**2 * (1e-4)**2 * spsfsontarget.val*spixelsperpsf.val * 4.66e6*1/1.12, "Cosmics"),
        #((1e-4/spixelsize.val)**2)*(4.66*1e6*1*1/1.12)*spsfsontarget.val*spixelsperpsf.val,"Cosmics"),
        ((10**sdarkcurrent.val)*spsfsontarget.val*spixelsperpsf.val,"DC"),
        (sreadnoise.val*spsfsontarget.val*spixelsperpsf.val,"RN")
        ],dt)
signalsax = plt.subplot(gs[17:21,10:15])
signalsax.set_ylabel(r"$N_{e^-}/sec$")
signalsax.set_ylim(0,1.1*signalsvals["flux"].max())
signalsrect = signalsax.bar(range(signalsvals.size),signalsvals["flux"],
                            align="center")
for r in signalsrect[1:]: r.set_color("red")
plt.xticks(np.indices(signalsvals["label"].shape)[0],
           signalsvals["label"])
# put the rectangles in a dictionary, so that you don't have to remember their order
signalsdict = dict((s,r) for s,r in zip(signalsvals["label"],signalsrect))

# SNR
time = np.linspace(1,3600,1000)
#snrax = plt.subplot(gs[14:19,9:16])
snrax = plt.subplot(gs[8:15,9:16])
snrax.set_xlim(0,time.max())
snrax.set_title("SNR")
snrax.set_xlabel("Integration time [min]")
snrax.set_xticks(np.linspace(0,3600,7))
snrax.set_xticklabels([int(i/60.) for i in snrax.get_xticks()])
snrax.ticklabel_format(style='sci',axis='y')
snrax.grid(True)
snrplt, = snrax.plot(time,
                     calcSNR(time=time, 
                             sig = signalsdict["Target"].get_height(),
                             bgnd = signalsdict["Sky"].get_height(),
                             Idc = signalsdict["DC"].get_height(),
                             readnoise = signalsdict["RN"].get_height()
                             )
                     )
ruinedax = snrax.twinx() # pixels ruined by cosmics
ruinedax.set_ylim(0,1)
ruinedax.set_xticks(np.linspace(0,3600,7))
ruinedax.set_xticklabels([int(i/60.) for i in ruinedax.get_xticks()])
ruinedax.set_ylabel("Pixel survival prob. from cosmics",labelpad=-1,rotation=-90,color='r')
singleprob, targetpixprob = calcCosmicHits(time=time)
ruinedplts = {}
ruinedplts["singlepix"], = ruinedax.plot(time,singleprob,'r--',label="single pixel")
ruinedplts["targetpix"], = ruinedax.plot(time,targetpixprob,'r-',label="target pixels")
legend = ruinedax.legend(loc=2,frameon=False)
for t in legend.get_texts():
    t.set_color("black")

# text formula
#fig.text(0.69,0.05,r"$\frac{S}{N}=\frac{St}{\sqrt{(S+Bn_{pix}+I_d n_{pix})t + R^2_n n_{pix}}}$",size="large")

## Dynamic updating
## signals bar chart
def update_signals_bar():
    signalsdict["Target"].set_height((10**stargets.val)*(np.pi*saperture.val**2)*sbandwidth.val*sqe.val)
    signalsdict["Sky"].set_height((10**ssky.val)*(np.pi*saperture.val**2)*sbandwidth.val*sqe.val)
    signalsdict["Cosmics"].set_height((10**scosmics.val) * (spixelsize.val)**2 * (1e-4)**2 * spsfsontarget.val*spixelsperpsf.val * 4.66e6*1/1.12),
    signalsdict["DC"].set_height((10**sdarkcurrent.val)*spsfsontarget.val*spixelsperpsf.val)
    signalsdict["RN"].set_height(sreadnoise.val*spsfsontarget.val*spixelsperpsf.val)
    # get maximum
    maxsig = max([rect.get_height() for rect in signalsdict.values()])
    signalsax.set_ylim(0,maxsig*1.1)
    #for r,s in zip(signalsrect,signalsvals["flux"]):
        #r.set_height(s)

def update_snr_plot():
    newsnr = calcSNR(time=time, 
                     sig = signalsdict["Target"].get_height(),
                     bgnd = signalsdict["Sky"].get_height(),    
                     Idc = signalsdict["DC"].get_height(),
                     readnoise = signalsdict["RN"].get_height(),
                     )
    maxsig = newsnr.max()
    snrax.set_ylim(0,maxsig*1.1)
    snrplt.set_ydata(newsnr)
    newcosmicsSingle, newcosmicsTarget=calcCosmicHits(time=time)
    ruinedplts["singlepix"].set_ydata(newcosmicsSingle)
    ruinedplts["targetpix"].set_ydata(newcosmicsTarget)

# whole figure
def update(val):
    ## signal sources
    mytargplt.set_xdata(stargets.val)
    myskyplt.set_xdata(ssky.val)
    mydarkcurrentplt.set_xdata(sdarkcurrent.val)
    mycosmicsplt.set_xdata(scosmics.val)
    myreadnoiseplt.set_xdata(sreadnoise.val)
    ## signals bar plot
    update_signals_bar()
    ## snr plot
    update_snr_plot()
    fig.canvas.draw_idle()

stargets.on_changed(update)
ssky.on_changed(update)
sdarkcurrent.on_changed(update)
scosmics.on_changed(update)
sreadnoise.on_changed(update)

saperture.on_changed(update)
sbandwidth.on_changed(update)
spixelsperpsf.on_changed(update)
spsfsontarget.on_changed(update)
spixelsize.on_changed(update)
sccdsize.on_changed(update)
sqe.on_changed(update)

# maximize window size
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())

plt.show()
