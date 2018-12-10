Title: check_current_status.sh

Purpose: check current status for running scanners such as which scanners are currently running, what image number the scanner is currently working and so on.

Instruction:

1. Open 'Terminal'

2. Set current directory that contains scripts
i.e. if 'Scripts' forlder contains scripts type line between " " in the command-line 
"cd /Desktop/Scripts"

3. Find out required information as follows
path: complete path that contains "*.log" files in other words the path that contains "scanner_api.sh" file
i.e. /home/spaldinglab/Desktop/production_Scripts/Scripts_main/

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

