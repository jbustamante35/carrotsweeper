#!/bin/bash

# $2: for first part type 1
#I don't know how to cutout the very first and last character from the string
#scanimage -L | grep Epson | cut -d' ' -f2
#`epkowa:interpreter:001:004'
# $1 is brand name for scanner connected 
# $3 device number

name=$(/usr/bin/time -f %e 2>&1 scanimage -L | grep $1 | grep $3 | cut -d'`' -f2 | cut -d':' -f$2)
#second=$(/usr/bin/time -f %e 2>&1 scanimage -L | grep $1 | cut -d'`' -f2 | cut -d':' -f$2)
echo $name
#scanner name part1
#scanimage -L | grep Epson | cut -d':' -f1
#scanner name part2
#scanimage -L | grep Epson | cut -d':' -f2

