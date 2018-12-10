#!/bin/bash
####################################
# run single scanner and save single image
# vars intput: 	$1: file name
#		$2: resolution
## Added April292015: variabilty for scanner manufacturer
#		$3: first part of scanner name 
#		$4: second part of scanner name 
#		$5: root hub
#		$6: device number
####################################
#echo -ne "$1 \n"
#echo -ne "\\************************** \n"
#echo -ne "CALL TO SINGLE SCAN CODE \n"
#echo -ne "file to save:" $1 "\n"
#echo -ne "scan resolution:" $2 "\n"
#echo -ne "usb root id:" $3 "\n"
#echo -ne "usb device number:" $4 "\n"
#echo -ne "************************** \n"
# run single scan with the inputs and record the time
## Original Code
#scanTime=$(/usr/bin/time -f %e 2>&1 scanimage --mode color -d "epson2:libusb:$3:$4" --resolution $2 --format tiff > $1)
## Made scanner name for scanimage code as variable to be passed
scanTime=$(/usr/bin/time -f %e 2>&1 scanimage --mode color -d "$3:$4:$5:$6" --x-resolution $2 --y-resolution $2 --format tiff > $1)
# compress image to new image	

# move compress image to old image

# delete compress old image


#echo -ne "************************** \n"
#echo "scan time:" $scanTime
echo $scanTime
#echo -ne "************************** \n"
