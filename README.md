### Planck

<p align="justify"> Planck is a rewrite of the previous project <a href="https://github.com/HemanthHaridas/plank.py"> plank.py</a>, which was a mixture of CPython and Cython. The older version suffered from poor performance due to its design and my limited programming skills at the time of its creation.</P>

<p align="justify"> Planck is instead written completely with a Python frontend (for easy installation and usage), with the computationally heavy routines implemented either purely using numpy, or in C++. Planck is not intended to be a replacement for commercially available codes like Gaussian, Orca and NWChem, but is a package that can be used for teaching introductory quantum chemistry for students.

<p align="justify"> One of the main goals of the Planck project is that it must be able to run from Jupyter Notebooks, which would be particulary handy for tutorial sessions.

