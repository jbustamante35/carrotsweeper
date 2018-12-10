#! /bin/bash
####################################
# update for 
# 1. Assume no other usb device is connected
# this method relies on two assumptions: 
# 1) that each line which contains the word vendor is a line which descibes a scanner
# 2) that the root hubs will always be listed with the second hub first and the first hub second

####################################
# LAPTOP specific version
# given the root hub and por tnumber return the device number
# vars intput: 	$1: root hub
#		$2: port number		
####################################
#!/bin/bash
# get the number of lines in the device tree

numports=$(lsusb -t |  grep -e "Port" -c)
numlines=$(($numports+2))
'
case $1 in

	001 )
		deviceID=$(lsusb -t | grep "Bus 01" -A$numlines | grep "Port 4" | cut -d ":" -f2 | cut -d "," -f1 | cut -d " " -f3);
		break;;
	002 ) 
		deviceID=$(lsusb -t | grep "Bus 02" -A$numlines | grep "Bus 01" -B$numlines | grep "Port 4" | cut -d ":" -f2 | cut -d "," -f1 | cut -d " " -f3);
		break;;
esac
'
if [ $1 == "002" ]; then
	deviceID=$(lsusb -t | grep "Bus 02" -A$numlines | grep "Bus 01" -B$numlines | grep "Port 4" | cut -d ":" -f2 | cut -d "," -f1 | cut -d " " -f3);
fi
# store the device ID of the scanner at the root hub $1 and the port $2
len=`echo $deviceID|wc -c`

# zero pad the decive number so that it is 3 charaters long
if [ $len -eq "2" ]; then
	deviceID="00"$deviceID;
fi

if [ $len -eq "3" ]; then
	deviceID="0"$deviceID;
fi
echo "$deviceID"




'
# search for vendor lines AFTER the key word Bus 01
if [ $1 == "001" ]; then
	deviceID=$(lsusb -t |  grep "Bus 01" -A$numlines | grep "Vendor" | grep "Port $2:" | cut -d "v" -f2 | cut -d "," -f1 | cut -d " " -f2);
fi
# search for vendor lines BEFORE the key word Bus 01
if [ $1 == "000" ]; then
	deviceID=$(lsusb -t |  grep "Bus 01" -B$numlines | grep "Vendor" | grep "Port $2:" | cut -d "v" -f2 | cut -d "," -f1 | cut -d " " -f2);
fi

# this line is for laptop
if [ $1 == "002" ]; then
	deviceID=$(lsusb -t | grep "Bus 02" -A$numlines | grep "Vendor" | grep "Port $2:" | cut -d "I" -f1 | cut -d ":" -f2 | cut -d "," -f1 | cut -d " " -f3);

fi

# this line is for laptop for multiple scanners
if [ $1 == "002" ]; then
	deviceID=$(lsusb -t | grep "Bus 02" -A20 | grep "Bus 01" -B20 | grep "Port 4" | cut -d ":" -f2 | cut -d "," -f1 | cut -d " " -f3);

fi
'



#Better way to grep device number?
#lsusb | grep Epson | cut -d ':' -f1 | cut -d ' ' -f4 





