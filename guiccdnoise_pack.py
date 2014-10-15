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

import warnings
import sys
import matplotlib as mpl
try:
    mpl.use('TkAgg')
except UserWarning:
    pass
if sys.version_info < (3,0):
    import Tkinter as Tk
    import Tkconstants
else:
    import tkinter as Tk
    import tkinter.constants as Tkconstants
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button#, RadioButtons
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import numpy as np


# some text formatting
mpl.rcParams["font.size"]=16
mpl.rcParams["text.color"]='white'
mpl.rcParams["axes.labelcolor"]='white'
mpl.rcParams["xtick.color"]='white'
mpl.rcParams["ytick.color"]='white'

# values for the sliders
xlims = (np.log10(1e-3), np.log10(1e10))
middle = np.mean(xlims)
targetcolor='blue'
bgndcolor='red'
telescopecolor='green'
sliderlength=400
sliderres=0.01 


class ccd_gui(Tk.Tk):
    def __init__(self, parent):
        Tk.Tk.__init__(self, parent)
        self.parent = parent
        self.initialize()
        self.initialize_target_options()
        self.initialize_background_options()
        self.initialize_telescope_options()
        self.initialize_signals_dict()
        self.initialize_snr_figure()
        self.initialize_signals_figure()
        self.initialize_layout()

    def initialize(self):
        #self.grid()
        #self.grid_columnconfigure(0,weight=1) # resize column 0 when window resized
        pass

    ### Sliders ###
    def initialize_target_options(self):
        self.starget=Tk.Scale(self,
                              from_=xlims[0],to=xlims[1], 
                              label='Targets [V-band mag]',
                              orient='horizontal',
                              troughcolor=targetcolor,
                              length=sliderlength,
                              resolution=sliderres,
                              showvalue=True,
                              tickinterval=(xlims[1]-xlims[0])/10,
                              command=self.OnSliderChange)

    def initialize_background_options(self):
        self.ssky = Tk.Scale(self, from_=xlims[0],to=xlims[1],
                             label='Sky',
                             orient='horizontal',
                             troughcolor=bgndcolor,
                             length=sliderlength,
                             resolution=sliderres)
        self.sdarkcurrent = Tk.Scale(self, from_=xlims[0],to=xlims[1], 
                                     label='Dark current',
                                     orient='horizontal',
                                     troughcolor=bgndcolor,
                                     length=sliderlength,
                                     resolution=sliderres)
        self.sreadnoise = Tk.Scale(self, from_=0, to=15,
                                   label='Read noise',
                                   orient='horizontal',
                                   troughcolor=bgndcolor,
                                   length=sliderlength,
                                   resolution=sliderres)
        self.scosmics = Tk.Scale(self, from_=-5, to=1,
                                 label='Cosmics',
                                 orient='horizontal',
                                 troughcolor=bgndcolor,
                                 length=sliderlength,
                                 resolution=sliderres)
                           
    def initialize_telescope_options(self):
        self.saperture=Tk.Scale(self, from_=0,to=10,
                                label= "Aperture diameter [m]",
                                orient='horizontal',
                                troughcolor=telescopecolor,
                                length=sliderlength,
                                resolution=sliderres)
        self.sbandwidth=Tk.Scale(self, from_=0,to=10,
                                 label= "Bandwidth [microns]",
                                 orient='horizontal',
                                 troughcolor=telescopecolor,
                                 length=sliderlength,
                                 resolution=sliderres)
        self.sccdsize=Tk.Scale(self, from_=500, to=5000,
                               label = "CCD dimension [pixels]",
                               orient='horizontal',
                               troughcolor=telescopecolor,
                               length=sliderlength,
                               resolution=sliderres)
        self.sqe=Tk.Scale(self, from_=0,to=1,
                          label="Quantum efficiency",
                          orient='horizontal',
                          troughcolor=telescopecolor,
                          length=sliderlength,
                          resolution=sliderres)
        self.spixelsize=Tk.Scale(self, from_=1,to=100,
                                 label="Pixel size [microns]",
                                 orient='horizontal',
                                 troughcolor=telescopecolor,
                                 length=sliderlength,
                                 resolution=sliderres)
        self.spixelsperpsf=Tk.Scale(self, from_=1,to=40,
                                    label="Pixels under PSF [#]",
                                    orient='horizontal',
                                    troughcolor=telescopecolor,
                                    length=sliderlength,
                                    resolution=sliderres)
        self.spsfsontarget=Tk.Scale(self, from_=1,to=200,
                                    label="Target size [PSFs]",
                                    orient='horizontal',
                                    troughcolor=telescopecolor,
                                    length=sliderlength,
                                    resolution=sliderres)
     
    def initialize_signals_dict(self):
        telescope_factor = np.pi*(self.saperture.get()**2)*\
            self.sbandwidth.get()*self.sqe.get()
        self.signals_dict = {'Target': (10**self.starget.get()) * telescope_factor ,
                             'Sky': (10**self.ssky.get()) * telescope_factor,
                             'Cosmics': (10**self.scosmics.get()) * 
                             (self.spixelsize.get())**2 * (1e-4)**2 * 
                             self.spsfsontarget.get()*self.spixelsperpsf.get() * 
                             4.66e6*1/1.12,
                             'DC': 10**self.sdarkcurrent.get() * 
                             self.spsfsontarget.get()*self.spixelsperpsf.get(),
                             'RN': self.sreadnoise.get() * 
                             self.spsfsontarget.get() *
                             self.spixelsperpsf.get()}

    def initialize_snr_figure(self):
        pass

    def initialize_signals_figure(self):
        self.f_signals = mpl.figure.Figure()
        self.ax_signals = self.f_signals.add_subplot(111)
        self.ax_signals.hist(np.random.normal(size=100))
        self.canvas_signals = FigureCanvasTkAgg(self.f_signals, master=self)

    def initialize_layout(self):
        # signal
        self.starget.pack(side=Tk.TOP, fill=Tk.X)
        # noise
        self.ssky.pack(side=Tk.TOP, fill=Tk.X)
        self.sdarkcurrent.pack(side=Tk.TOP, fill=Tk.X)
        self.sreadnoise.pack(side=Tk.TOP, fill=Tk.X)
        self.scosmics.pack(side=Tk.TOP, fill=Tk.X)
        # telescope
        self.saperture.pack(side=Tk.TOP, fill=Tk.X)
        self.sbandwidth.pack(side=Tk.TOP, fill=Tk.X)
        self.sccdsize.pack(side=Tk.TOP, fill=Tk.X)
        self.sqe.pack(side=Tk.TOP, fill=Tk.X)
        self.spixelsize.pack(side=Tk.TOP, fill=Tk.X)
        self.spixelsperpsf.pack(side=Tk.TOP, fill=Tk.X)
        self.spsfsontarget.pack(side=Tk.TOP, fill=Tk.X)
        # figures
        self.canvas_signals.get_tk_widget().pack(side=Tk.RIGHT)
        

    ### Actions
    def OnSliderChange(self, value):
        """
        Update the figures using the slider values
        """
        pass

if __name__ == "__main__":
    gui = ccd_gui(None)
    gui.title('CCD Noise!')
    gui.mainloop() # start looping, looking for user interaction
