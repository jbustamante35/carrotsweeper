#!/bin/bash
####################################
# run single scanner and save single image
# vars intput: 	$1: file name
#		$2: resolution
#		$3: root hub
#		$4: device number
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
scanTime=$(/usr/bin/time -f %e 2>&1 scanimage --mode color -d "epson2:libusb:$3:$4" --resolution $2 --format tiff > $1)
# compress image to new image	

# move compress image to old image

# delete compress old image


#echo -ne "************************** \n"
#echo "scan time:" $scanTime
echo $scanTime
#echo -ne "************************** \n"
