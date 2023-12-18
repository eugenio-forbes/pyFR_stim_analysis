%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Behavioral Toolbox Release 1.1
          July 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Table of Contents:
- Introduction
- Analyses
  - Input format
  - Output format
- Helpers
  - Masks
  - matrixops
  - structs
  - transitions
- Plot
- Contact

============
INTRODUCTION
============

Welcome to the Behavioral Toolbox Release 1! 
This first release has functions designed for performing basic analyses on free recall studies. This file is meant to orient you to the main structure of these functions: what they expect as input and what they produce as output. Although much of this information is contained in each of the help descriptions for each function, you may find it useful to understand the common themes of the functions by reading through this document.
This file is organized according to the directories of functions that comprise the Behavioral Toolbox Release 1.
These scripts were written for use with Matlab version 7.5 (R2007b), but should work with any recent version of Matlab. The directory and subdirectories containing the behavioral toolbox functions must be in your Matlab path in order for these functions to work properly.

========
ANALYSES
========

Most of these functions expect the associated information from a free recall study to be in a particular format. We give an overview of this format here.

*** Inputs ***
In a free recall study, each trial has certain information associated with it: the items presented, the items recalled, the corresponding subject.
Combined across all such trials, one can generate matrices where each row represents a trial, each column represents a recalled item. Specifically, column i represents output position i, and if no item was recalled for that output position, it is simply left as 0. Two main ways of representing these recalls are to index recalled items by the serial position for that presented list: integers from 1 to the list-length (referred to as *recalls_matrix* in many toolbox functions). Another way is to index items by their number in the word pool (referred to as *rec_itemnos*). This allows for more detailed information to be extracted, such as if any items recalled were prior-list intrusions.
Corresponding to each trial row, one can also generate matrices with each of the items presented according to their number in the wordpool (pres_itemnos)
Critical to most functions is a vector where each row corresponds to the number of the subject for that trial (referred to as *subjects* in toolbox functions).
To facilitate analyses of extra-list and prior-list intrusions, many functions expect as input a matrix with intrusion information (called *intrusions*). For the specifics of how this matrix can be created and is designed, see make_intrusions.m.

*** Outputs ***

All of the functions in this release output a number or a set of numbers for each participant, depending on the analysis. When the analysis outputs one number per participant (e.g. proportion of items recalled), the output is a vector, where each row corresponds to the number for one participant. The numbers are listed according to the ascending order of subject numbers. When the analysis outputs more than one number per subject, the rows still correspond to the ascending order of subject numbers, and the columns index the different values for the analysis for that one subject (e.g. in spc, the columns indicate probability of recall at ascending serialpositions).

=======
HELPERS
=======

Many of the functions in the Helpers subdirectories are internal functions meant to supplement the functions provided in the Analyses folder. For Release 1, only the functions in the Masks directory and the function make_intrusions.m in the Matrixops folder are meant are designed to be explicitly called upon by the user. make_intrusions.m was described above in the Inputs section of this file, and Masks are described in more detail below.
If you would like to use any of the other functions on their own, the detailed docstrings should suffice. Good luck!

*** Masks ***

Masks are matrices with logical elements used to 'mask' out particular recalled and/or presented items. For instance, suppose we are only interested in recall of the odd-numbered items. One could simply create a mask the same size as pres_itemnos (see Input Format) where all even-numbered items are false, and all odd-numbered items are true. One could then create a similar mask for the recalled items, and then pass these two masks into the appropriate analysis function.
NOTE: Many functions make default masks if they are not provided by the user.

====
PLOT
====

The plot functions have been designed especially to take the output of functions from the Analyses folder of the Behavioral Toolbox Release 1, and turn them into lovely plots. Simply pass in the numbers from those functions, and admire the figures! More detailed information about how to use those functions are provided in their docstrings.

- plot_crp.m can plot the output of lag_crp.m
- plot_spc.m can plot the output of spc.m OR pfr.m 

=======
CONTACT
=======

Behavioral Toolbox Release 1 was created at the Computational Memory Laboratory of Dr. Michael Kahana at the University of Pennsylvania. For questions concerning the questions contained here, please contact ...
