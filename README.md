# xraytools
Simple algorithms for loading and processing small-angle X-ray scattering
and X-ray diffraction data files.
The two main resource here are:
- reader functions for XRDML (like Malvern Panalytical RXD devices) and Xenocs data files
- tools to plot, subtract data, do background correction and peak finding

The smoothing, background correction and peak finding are part of a more general signal
package 'tomio.signals'.

The plotting funcitons are not meant as an extension of the plot() methode, because they focus
on the specific type of x-ray data. It should not break anything, but it may raise a warning from
the build tools.
