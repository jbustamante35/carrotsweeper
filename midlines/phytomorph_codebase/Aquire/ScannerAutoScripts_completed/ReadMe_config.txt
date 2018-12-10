Title: config.txt

Purpose: Set up a config file for running scanners such as resolution, total duration and so on.

Instruction:
User input configuration 
intput: 	resolution: X in dpi, desired resolution
		example: resolution:1200

		deltaTime: X in seconds, desired break time for each scan
		example: deltaTime:300

		duration: X in seconds, desired total time for running scanner
		example: duration:86400

		imagepath: complete path to store the images that are scanned
		example: imagepath:/home/spaldinglab/Desktop/images

		scannerID: (RootHub in 3 digits)_(port number in a single digit)
		scannerID example: root hub-2, port-4 ==> 002_4
scannerID4:002_4


scannerID4:002_4
scannerID1:001_3
scannerID5:002_5
scannerID6:002_6
scannerID2:002_2
scannerID3:002_3
1. Open 'Terminal'

2. Set current directory that contains scripts
i.e. if 'Scripts' forlder contains scripts type line between " " in the command-line 
"cd /Desktop/Scripts"

3. Find out required information as follows
path: complete path that contains "*.log" files in other words the path that contains "scanner_api.sh" file
i.e. /home/spaldinglab/Desktop/Scripts/

P.S. "*.log" includes (1.status, 2. log) for each scanner run

4. Type following line in the command-line to run scanner 
Replace [path] specific information from 3
"./check_current_status.sh [path]"

5. Repeat 4 to check current status of working scanners

ReadMe_written_Date: March.24.2015.
Code_Modified_Date:  
Version:
Author:
Contact:
License: 

