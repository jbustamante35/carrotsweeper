#! /bin/bash
####################################
# given the root hub and por tnumber return the device number
# vars intput: 	$1: root hub
#		$2: port number		
####################################
# get the number of lines in the device tree
numlines=$(lsusb -t |  grep -e "Port" -c)
# this method relies on two assumptions: 
# 1) that each line which contains the word vendor is a line which descibes a scanner
# 2) that the root hubs will always be listed with the second hub first and the first hub second

# search for vendor lines AFTER the key word Bus 01
if [ $1 == "001" ]; then
	deviceID=$(lsusb -t |  grep "Bus 01" -A$numlines | grep "Vendor" | grep "Port $2:" | cut -d "v" -f2 | cut -d "," -f1 | cut -d " " -f2);
fi
# search for vendor lines BEFORE the key word Bus 01
if [ $1 == "002" ]; then
	deviceID=$(lsusb -t |  grep "Bus 01" -B$numlines | grep "Vendor" | grep "Port $2:" | cut -d "v" -f2 | cut -d "," -f1 | cut -d " " -f2);
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


#Better way to grep device number?
#lsusb | grep Epson | cut -d ':' -f1 | cut -d ' ' -f4 





