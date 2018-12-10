#!/bin/bash
# April 27 2015
# This script will help to create configuration txt file for running scanner.
## 1) Add this to head of .m file #!/usr/bin/octave -qf 
## 2) Print .m file in command line $ vi .m
## 3) Make .m file executable in command line $ chmod 755 .m  
## 4) Run $ ./ .m
## 5) Do this for all related .m files

#/home/spaldings_lab/Desktop/Scripts/csv2struct.m
#/home/spaldings_lab/Desktop/Scripts/readtext.m
#/home/spaldings_lab/Desktop/Scripts/selectComputerFromList.m
#/home/spaldings_lab/Desktop/Scripts/selectFileFromTitle.m
# This line works
# Run .m code
cd /home/spaldings_lab/Desktop/Scripts 
#echo pwd
octave CL_interface.m 
#| save -ascii test35.txt
#1>test34.txt
# Sole parts of CL_interface.m works but it does not call other *.m scripts related.
# I chmod for all *.m but still it does not work.
# Variables work itself? or store standard out to txt?
## April 28 2015
# 1>standardout.txt or save standardout.txt do not seem to work here 
