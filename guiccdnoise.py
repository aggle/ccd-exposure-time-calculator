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
- Run from ipython with %run guiccdnoise.py or from command line with ./guiccdnoise.py
- Click on the sliders to set the parameter values. 


to-do
----------
- Make the display as useful as the matplotlib version
"""

import matplotlib as mpl
try:
    mpl.use('TkAgg')
except UserWarning:
    pass

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button#, RadioButtons
import warnings
import sys

if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg


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

dt = np.dtype([("flux",'f8'),("label","S20")])

# units: Photons/(um m^2 sec)
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

# Sky, Photons/(um m^2 sec)
sky = np.array([
        (6.1e3,"Full moon"), 
#        (5.0e2, "Half moon"), 
        (2.9e2, "New moon"),
        ], dt)
sky.sort(order="flux")

# Dark current, e-/(pixel sec)
darkcurrent = np.array([
        (6.2e-3, "HST WFC"),
        (8e-1, "Keck NIRSPEC"),
        (1e-1, "Commercial CCD")
#        (1e-3,"40 C"), 
#        (1e-1,"80 C"), 
#        (1e2, "100 C")
        ], dt)
darkcurrent.sort(order="flux")

# Read noise, e-/pix/read
readnoise = np.array([
        (2., "Chandra ACIS"),
        (3., "HST WFC3"),
        (4., "P1640"),
        (10., "Subaru Suprime-Cam")
        ], dt)
readnoise.sort(order="flux")

# Cosmics, N/(cm^2 sec)
cosmics = np.array([
        (1.e-2,"Sea lvl"),
        (2e-2, "4 km"),
        (1,"Space"),
        (1.3e-4, "LHC")
        ], dt)
cosmics.sort(order="flux")

# some slider values
xlims = (np.log10(1e-3), np.log10(1e10))
middle = np.mean(xlims)
targetcolor='blue'
bgndcolor='red'
telescopecolor='green'
sliderlength=400
sliderres=0.01 

root = Tk.Tk()
root.wm_title("CCD SNR calculator")

f_snr = mpl.figure.Figure()#figsize=(5,4), dpi=100)
f_signals = mpl.figure.Figure()#figsize=(5,4), dpi=100)
canvas_snr = FigureCanvasTkAgg(f_snr, master=root)
canvas_signals = FigureCanvasTkAgg(f_signals, master=root)
canvas_snr.show()
canvas_signals.show()

def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate

### Sliders
starget = Tk.Scale(root, from_=xlims[0],to=xlims[1], 
                   label='Targets [V-band mag]',
                   orient='horizontal',
                   troughcolor=targetcolor,
                   length=sliderlength,
                   resolution=sliderres,
                   tickinterval=1)
starget.set(np.log10(targets['flux'][0]))
ssky = Tk.Scale(root, from_=xlims[0],to=xlims[1],
                     label='Sky',
                     orient='horizontal',
                     troughcolor=bgndcolor,
                     length=sliderlength,
                     resolution=sliderres)
ssky.set(np.log10(sky['flux'][0]))
sdarkcurrent = Tk.Scale(root, from_=xlims[0],to=xlims[1], 
                             label='Dark current',
                             orient='horizontal',
                             troughcolor=bgndcolor,
                             length=sliderlength,
                             resolution=sliderres)
sdarkcurrent.set(np.log10(darkcurrent['flux'][0]))
sreadnoise = Tk.Scale(root, from_=0, to=15,
                           label='Read noise',
                           orient='horizontal',
                           troughcolor=bgndcolor,
                           length=sliderlength,
                           resolution=sliderres)
sreadnoise.set(np.log10(readnoise['flux'][0]))
scosmics = Tk.Scale(root, from_=-5, to=1,
                         label='Cosmics',
                         orient='horizontal',
                         troughcolor=bgndcolor,
                         length=sliderlength,
                         resolution=sliderres)
scosmics.set(np.log10(cosmics['flux'][0]))
### Telescope parameters
saperture=Tk.Scale(root, from_=0,to=10,
                        label= "Aperture diameter [m]",
                        orient='horizontal',
                        troughcolor=telescopecolor,
                        length=sliderlength,
                        resolution=sliderres)
saperture.set(1)
sbandwidth=Tk.Scale(root, from_=0,to=10,
                         label= "Bandwidth [microns]",
                         orient='horizontal',
                         troughcolor=telescopecolor,
                         length=sliderlength,
                         resolution=sliderres)
sbandwidth.set(1)
sccdsize=Tk.Scale(root, from_=500, to=5000,
                       label = "CCD dimension [pixels]",
                       orient='horizontal',
                       troughcolor=telescopecolor,
                       length=sliderlength,
                       resolution=sliderres)
sccdsize.set(2048)
sqe=Tk.Scale(root, from_=0,to=1,
                  label="Quantum efficiency",
                  orient='horizontal',
                  troughcolor=telescopecolor,
                  length=sliderlength,
                  resolution=sliderres)
sqe.set(0.85)
spixelsize=Tk.Scale(root, from_=1,to=100,
                         label="Pixel size [microns]",
                         orient='horizontal',
                         troughcolor=telescopecolor,
                         length=sliderlength,
                         resolution=sliderres)
spixelsize.set(18)
spixelsperpsf=Tk.Scale(root, from_=1,to=40,
                            label="Pixels under PSF [#]",
                            orient='horizontal',
                            troughcolor=telescopecolor,
                            length=sliderlength,
                            resolution=sliderres)
spixelsperpsf.set(10)
spsfsontarget=Tk.Scale(root, from_=1,to=200,
                            label="Target size [PSFs]",
                            orient='horizontal',
                            troughcolor=telescopecolor,
                            length=sliderlength,
                            resolution=sliderres)
spsfsontarget.set(1)

### Signal-to-Noise
def calcSNR(sig=1,bgnd=0,readnoise=0,Idc=0,npix=1,time=1):
    """
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

### Cosmics
def calcCosmicHits(time=1):
    """
    The probability of losing a pixel to cosmics
    """
    time = np.asarray(time)
    ccdflux = (spixelsize.get()*sccdsize.get()*1e-4)**2 * 10**scosmics.get() # flux through ccd
    npix = spsfsontarget.get()*spixelsperpsf.get()
    probPerPix = (1-1./(sccdsize.get()**2))**(ccdflux*time)
    probPerTargetPix = probPerPix**npix
    return (probPerPix,probPerTargetPix)

### Signals bar chart
signalsvals = np.array([
        (mag2flux(starget.get())*(np.pi*saperture.get()**2)*sbandwidth.get()*sqe.get(),"Target"),
        (mag2flux(ssky.get())*(np.pi*saperture.get()**2)*sbandwidth.get()*sqe.get(),"Sky"),
        (mag2flux(scosmics.get()) * (spixelsize.get())**2 * (1e-4)**2 * spsfsontarget.get()*spixelsperpsf.get() * 4.66e6*1/1.12, "Cosmics"),
        (mag2flux(sdarkcurrent.get())*spsfsontarget.get()*spixelsperpsf.get(),"DC"),
        (sreadnoise.get()*spsfsontarget.get()*spixelsperpsf.get(),"RN")
        ],dt)
signalsax = f_signals.add_subplot(111)
signalsax.set_ylabel(r"$N_{e^-}/sec$")
signalsax.set_ylim(0,1.1*signalsvals["flux"].max())
signalsax.set_xticks(np.indices(signalsvals["label"].shape)[0])
signalsax.set_xticklabels(signalsvals["label"])
signalsrect = signalsax.bar(range(signalsvals.size),signalsvals["flux"],
                            align="center")
signalsrect[0].set_color("blue")
for r in signalsrect[1:]: r.set_color("red")
# put the rectangles in a dictionary, so that you don't have to remember their order
signalsdict = dict((s,r) for s,r in zip(signalsvals["label"],signalsrect))

### SNR plot
time = np.linspace(1,3600,1000)
snrax = f_snr.add_subplot(111)
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
ruinedax = snrax.twinx()
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

def update_signals_bar():
    signalsdict["Target"].set_height(mag2flux(starget.get())*(np.pi*saperture.get()**2)*sbandwidth.get()*sqe.get())
    signalsdict["Sky"].set_height(mag2flux(ssky.get())*(np.pi*saperture.get()**2)*sbandwidth.get()*sqe.get())
    signalsdict["Cosmics"].set_height(mag2flux(scosmics.get()) * (spixelsize.get())**2 * (1e-4)**2 * spsfsontarget.get()*spixelsperpsf.get() * 4.66e6*1/1.12),
    signalsdict["DC"].set_height(mag2flux(sdarkcurrent.get())*spsfsontarget.get()*spixelsperpsf.get())
    signalsdict["RN"].set_height(sreadnoise.get()*spsfsontarget.get()*spixelsperpsf.get())
    # get maximum
    maxsig = max([rect.get_height() for rect in signalsdict.values()])
    signalsax.set_ylim(0,maxsig*1.1)
    #for r,s in zip(signalsrect,signalsvals["flux"]):
    #    r.set_height(s)

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
    ## signals bar plot
    update_signals_bar()
    ## snr plot
    update_snr_plot()
    f_snr.canvas.draw_idle()
    f_signals.canvas.draw_idle()

# add update() to the sliders
starget.config(command=update)
ssky.config(command=update)
sdarkcurrent.config(command=update)
sreadnoise.config(command=update)
scosmics.config(command=update)
### Telescope parameters
saperture.config(command=update)
sbandwidth.config(command=update)
sccdsize.config(command=update)
sqe.config(command=update)
spixelsize.config(command=update)
spixelsperpsf.config(command=update)
spsfsontarget.config(command=update)


### Placement
# noise sliders
noisecol=1
starget.grid(row=1,column=noisecol)
ssky.grid(row=2,column=noisecol)
sdarkcurrent.grid(row=3,column=noisecol)
sreadnoise.grid(row=4,column=noisecol)
scosmics.grid(row=5,column=noisecol)
# telescope sliders
telcol=1
saperture.grid(column=telcol)
sbandwidth.grid(column=telcol)
sccdsize.grid(column=telcol)
sqe.grid(column=telcol)
spixelsize.grid(column=telcol)
spixelsperpsf.grid(column=telcol)
spsfsontarget.grid(column=telcol)
# matplotlib figures
graphcol=3
canvas_snr.get_tk_widget().grid(row=1, column=graphcol, rowspan=6,columnspan=3,
                                sticky='W'+'E'+'N'+'S')
canvas_signals.get_tk_widget().grid(row=7, column=graphcol, rowspan=12,columnspan=3,
                                    ipadx=0,ipady=0)

if __name__ == "__main__":
    try:
        Tk.mainloop()
    except KeyboardInterrupt:
        sys.exit()
    #root.mainloop()
