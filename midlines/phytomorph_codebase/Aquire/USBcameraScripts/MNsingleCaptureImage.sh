#!/bin/bash
pkill gvfs
D=$(date +"%m-%d-%y")
folder=$HOME'/maizeData/seedlingData/'$D'/'
echo $folder
mkdir -p $folder
tmpLocation=$HOME'/tmpData/'
mkdir -p $tmpLocation
gphoto2 --capture-image-and-download --force-overwrite --filename=$tmpLocation/tmpImage.nef

