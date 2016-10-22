This is a suite of EBSD analysis software for 2D data sets collected
on a Quanta machine using Oxford software. While this is by no means
the fastest code of its kind, I wrote it so I would have a pared down
way to analyze EBSD data for my specific data sets. This is a pretty 
short set of code that includes some basic maps and pole figures and
statistical analysis using Kent distributions.

The interface is in master interface.nb. Different features are 
embedded in different functions to help streamline the process, but
converting source files to .meap files will be by far the longest 
process. This code is a work in progress in terms of making it general
enough for others to use; some of the exporting and importing features
are designed around my file nomenclature scheme, and the code is 
designed for aragonite data sets with 2 ROI.

Some auxiliary functions are in misorientation and kent_sp, but the
code that controls the interfaces is in masterFunctions.