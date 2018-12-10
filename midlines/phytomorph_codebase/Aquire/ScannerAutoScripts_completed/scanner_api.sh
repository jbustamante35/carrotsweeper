#!/bin/bash
####################################
# User input configuration 
# intput: 	$1: scannerID
#		$2: config.txt
#		$3: $this_computer from Run.sh
#		duration: X in seconds
#		imagepath: example /home/spaldinglab/Desktop/images
#		scannerID: (RootHub in 3 digits)_(port number in a single digit)
#		scannerID example: root hub-2, port-4 ==> 002_4
## April.8.2015.
### April.14.2015.
####################################
####################################
## Look for device number to pass on to run_scanner.sh
## how to feed text file in shell script??
## path that to look for ScannerWiring.txt. this_computer variable contains computer number and it was 
## found from Run.sh

path=/mnt/snapper/kernelSwellingData/config/$3
echo $path

## store ScannerWiring.txt into line

# init scann count
scannerCount=1
# This will read lines as long as we have some left
### I need to input data at the end since 1) I don't know how to deal with 2 txt in put 2) If I do 
### cat $path/ScannerWiring.txt | while read line, then somehow the variable is restored after the loop
while read line          
do           
	# read a scanner line
    	temp=$(echo $line | cut -d: -f2)
    	pretemp=$(echo $line | cut -d: -f1)

	# get the true scanner ID
	truescannerID[$scannerCount]=$(echo $pretemp | cut -d"D" -f2 | cut -d":" -f1)

	# get the scanner ID
	scannerID[$scannerCount]=$(echo $pretemp | cut -d"D" -f2)

	# get the scanner root
	scannerIDroot[$scannerCount]=$(echo $temp | cut -d_ -f1)

	# get the scanner port
	scannerIDport[$scannerCount]=$(echo $temp | cut -d_ -f2)
	# query for the scanner device number
	scannerIDdevnum[$scannerCount]=$(/home/spaldings_lab/Desktop/Scripts/Scripts_main/getDeviceForPort.sh ${scannerIDroot[$scannerCount]} ${scannerIDport[$scannerCount]})
	# report the found scanner id, root usb card, usb port, device number	
	echo "Request scanner ID:" ${scannerID[$scannerCount]}  " at root:" ${scannerIDroot[$scannerCount]} ":port" ${scannerIDport[$scannerCount]} ":device number" ${scannerIDdevnum[$scannerCount]}
	# increase the count to read into array
	((scannerCount++))
done <$path/ScannerWiring.txt

for index in `seq 1 ${#scannerIDroot[@]}`;
do
if [ $1 == ${truescannerID[$index]} ]; then
	j=$index
fi
done

echo "Scanning at scanner number" $1 ${scannerIDroot[$j]} ":" ${scannerIDport[$j]}

#### it seems that lost variable after while loop
## scannerIDdevnum contains device number. Find a way to match this info with input scanner number and output 
## single device number to run a scanner.
#### How to feed variable to 'run_scanner.sh' without changing class/structure of the variable.
#### How to deal with scannerCount. It always add up 1 more than what it should be at the end


#echo $scannerCount
#printf ${scannerIDdevnum[@]}
#printf ${truescannerID[@]}
#printf ${scannerID[@]}
#printf ${scannerIDroot[@]}
#printf ${scannerIDport[@]}
dev=$(echo ${scannerIDdevnum[$j]})
truenum=$(echo ${truescannerID[$j]})
ID=$(echo ${scannerID[$j]})
root=$(echo ${scannerIDroot[$j]})
port=$(echo ${scannerIDport[$j]})
#dev= ${scannerIDdevnum[$j]}
#truenum=${truescannerID[$j]}
#ID=${scannerID[$j]}
#root=${scannerIDroot[$j]}
#port=${scannerIDport[$j]}

##############################################
#changes 'date' format by replacing ' ' to '_' and ':' to '-'
#1)It is not good to have ' ' in the file/folder name
#2)snapper does not accept ':' in the file/folder name
##status_name=$(echo "$(date)" | sed -e 's/ /_/g'| sed -e 's/:/-/g')
#name for status text file
##status_name="status_file_"$status_name
#name for log file
##stdout_name="stdout_file_"$status_name
#takes logs into text file? and call 'run_scanner.sh'
/home/spaldings_lab/Desktop/Scripts/Scripts_main/run_scanner.sh $1 $status_name < $2 1 > $stdout_name.txt $dev $truenum $ID $root $port




