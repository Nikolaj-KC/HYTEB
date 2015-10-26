Welcome to HYTEB's documentation!
=================================
This is the official documentation for the ``HYTEB``.
This document aims at guiding the documenting and the usage.

What is HYTEB?
--------------
``HYTEB`` is a set of python modules for user-friendly software for analysing the data worth of geophysics (electromagnetic data) for improving the predictive capability of groundwater models.

Authors
-------
``HYTEB`` is developed by Nikolaj Kruse Christensen as part as his Ph.D. project, under
supervision of Steen Christensn, all of the Department
of Geoscience, Aarhus University, Denmark. The author welcomes interested third
party developers. This document is a work in progress.

Contact
-------
Nikolaj Kruse Christensen - http://pure.au.dk/portal/da/persons/id(9402400e-29ed-465f-904d-85c0c49e217e).html/ 

mail to - nikolaj.kruse@geo.au.dk

Documentation
-------------
The goal is selv-documenting code in sphinx, so no manual will be available.
The documentation is by default available in the html.


Get started with HYTEB
----------------------
We have only uploaded software that we developed. Software developed by others, for example AarhusInv, PEST, MODFLOW, TPROGS, BLOCKSIS, etc. must be purchased or downloaded from websites of the respective developers.

HYTEB is written in the python 2.7 syntax.
I recommend the Anaconda scientific python distribution , which includes the dependencies for ``HYTEB``. 

[https://store.continuum.io/cshop/anaconda/](https://store.continuum.io/cshop/anaconda/])

Once installed, clone (or download) the ``HYTEB`` repository. 

Examples
--------

See the paper in HESSD for a description of the elements in HYTEB and examples for:

http://www.hydrol-earth-syst-sci-discuss.net/12/9599/2015/hessd-12-9599-2015.pdf

In the paper 3 different inversion approaches are compared (Hydrogeological inversion, sequential hydrogeophysical inversion, joint hydrogeophysical inversion). 

For an example for the HYTEB setup for doinf joint hydrogeophysical inversion see: [main_HYTEB_joint_inversion.py](https://github.com/Nikolaj-KC/HYTEB/blob/master/main_HYTEB_joint_inversion.py) 


Links
=====
[AarhusInv - http://hgg.au.dk/software/aarhusinv/](http://hgg.au.dk/software/aarhusinv/)

[BlockSIS - http://www.iamg.org/documents/oldftp/VOL32/v32-10-12.zip/](http://www.iamg.org/documents/oldftp/VOL32/v32-10-12.zip)

[MODFLOW - http://water.usgs.gov/ogw/modflow/MODFLOW.html](http://water.usgs.gov/ogw/modflow/MODFLOW.html)

[PEST - http://www.pesthomepage.org/](http://www.pesthomepage.org/)

[T-PROGS - http://www.aquaveo.com/software/gms-tprogs](http://www.aquaveo.com/software/gms-tprogs)

References
==========
[Auken et al. 2014 - An overview of a highly versatile forward and stable inverse algorithm for airborne, ground-based and borehole electromagnetic and electric data](http://hgg.au.dk/fileadmin/www.gfs.au.dk/AUKEN2014D.pdf)

[Auken et al. 2015 - Manual for the inversion program AarhusInv](http://www.hgg.geo.au.dk/HGGsoftware/em1dinv/em1dinv_manual.pdf)

Doherty, J., 2010a, PEST, Model-independent parameter estimation—User manual (5th ed., with slight additions):
Brisbane, Australia, Watermark Numerical Computing.

Doherty, J., 2010b, Addendum to the PEST manual: Brisbane, Australia, Watermark Numerical Computing.



Publications
============
``HYTEB`` has been used to produce results in the following scientific
publications and presentations:

Christensen, N. K., Christensen, S., and Ferre, T. P. A.: A framework for testing the use of electric and electromagnetic data to reduce the prediction error of groundwater models, Hydrol. Earth Syst. Sci. Discuss., 12, 9599-9653, doi:10.5194/hessd-12-9599-2015, 2015
