iHSComputer
===========

computes iHS statistic accounting for missing data.

iHSComputer is written by Sandeep Venkataram and Yuan Zhu, Petrov Lab, Stanford University Department of Biology
Two versions of the software are provided. 
The default program is single threaded and requires gcc with c++98. The multithreaded version requires gcc with c++0x and pthread support.


To compile the single threaded version, just execute the command "make iHS" in the directory this archive was extracted to. 
To compile the multithreaded version, execute the command "make iHS_multithreaded". 
To compile the normalization program, execute the command "make normalizer"
Precompiled linux binaries are provided for your use