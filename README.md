ccd-exposure-time-calculator
============================

Usage:
=======
- Clone into the repository of your choice
- From a bash terminal, run: python guiccdnoise.py
- From an ipython termina, run: %run guiccdnoise.py

A window will appear with a number of sliders on the right and two 
Set the sliders to match the properties of your observations. The plots will rescale themselves to reflect how your choices affect the exposure time needed to overcome different sources of noise.

Sliders
=======
- V-band flux [mag]: V-band magnitude of your target
- Sky flux [mag]: V-band magnitude of the sky
- Dark current [units]: Units?
- Read noise [e/s]: read noise per pixel
- Cosmics [units]: units? density of cosmic rays?
- Aperture diameter [m]: mirror size
- Bandwidth [um]: effective filter width
- CCD dimension [pixels]: 
- Quantum efficiency [0-1]
- Pixel size
- Pixels under PSF [# pixels]
- Target size [# PSFs]

Plots
=====
Top: 
- the LHS shows the SNR as a function of time based on the noise, target, and instrument properties set on the sliders
- the RHS shows the probability to have a cosmic impact a single pixel based on the cosmic flux density

Bottom: histogram electrons per second broken down by source
