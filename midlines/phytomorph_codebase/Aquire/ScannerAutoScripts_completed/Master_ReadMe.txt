Title: Complete process flow

Flow: 
Scanning
A) 1. Scott_run.sh --> 2. (scanner_api.sh, create stdout*.txt) --> 3. run_scanner.sh --> 4. (getDeviceForPort.sh, run_single_scanner.sh, create *.log) --> 5. scanimage

Checking 
B) check_current_status.sh

Terminating
C) kill_all_scanner.sh

Idea: 
A)
1. 
User run 'Scott_run.sh' with the scanner number that is written on each scanner.

--'Scott_run.sh' takes scanner ID. 'config.txt' is hard coded

2. 'Scott_run.sh' calls a) 'scanner_api.sh' in 'production_scripts/Scripts_main/' with given scanner number in 1 and specified 'config.txt' in the same directory and creates b) 'stdout*.txt' for 'check_current_status.sh'.

--'scanner_api.sh' takes scanner ID and 'config.txt' from 'Scott_run.sh' command. Then it calls 'run_scanner.sh' and creates 'stdout*.txt'.

3. 'scanner_api.sh' calls 'run_scanner.sh', the core code to run scanner with given information in the previous steps.

--'run_scanner.sh' includes a)creating folder for images b)getting necessary information from the computer-device number and such c)compressing images and deleting the original images d)creating '*.log'.

4. 'run_scanner.sh' calls a) 'getDeviceForPort.sh' to get device number for the scanner which changes every on/off and 
b) 'run_single_scanner.sh' to run a single scanner and creates c) '*.log' for 'check_current_status.sh'.

--a)'getDeviceForPort.sh' search connected devices and get proper lines that includes device numbers for the scanners. Then add proper number of 0's for device ID. b)'run_single_scanner.sh' runs one scanner with given information from previous stages. "epson2:libusb:$3:$4" specifies the scanner that are using. We may change it if we replace current scanners to something else. Scanning time is measured in this script. c) '*.log' is created and stores in the same folder that contains scripts.

5. 'run_single_scanner.sh' calls command 'scanimage' with given configuration in the previous statges.

B)
Read log from '*.log' which stores all the log created for each scanner. Find out correct scanner by using starting time of the scanner. This staring time varies by each and every scanner.

C)
Kill all the processes related to the scanner code.
It works with name of the process

ReadMe_written_Date: March.26.2015.
Code_Modified_Date1: April.6.2015.
Code_Modified_Date2: 
Version:
Author:
Contact:
License:


#Information on Geometry for scanimage
  Geometry:
    -l 0..215.9mm [0]
        Top-left x position of scan area.
    -t 0..297.18mm [0]
        Top-left y position of scan area.
    -x 0..215.9mm [215.9]
        Width of scan-area.
    -y 0..297.18mm [297.18]
        Height of scan-area.
