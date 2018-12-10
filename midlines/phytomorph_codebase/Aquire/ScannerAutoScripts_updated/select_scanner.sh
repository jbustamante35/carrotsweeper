#!/bin/bash
# this is for selecting scanner and running scanner in a single terminal
echo "Which scanner do you want to run or EXIT?"
path2config=
path2run=
select selection in "Scanner1" "Scanner2" "Scanner3" "Scanner4" "Scanner5" "Scanner6" "Change_Option" "EXIT"; do
	case $selection in
		Scanner1 )
			SD=1;
			break;;
		Scanner2 ) 
			SD=2;
			break;;
		Scanner3 )
			SD=3;
			break;;
		Scanner4 ) 
			SD=4;
			break;;
		Scanner5 )
			SD=5;
			break;;
		Scanner6 ) 
			SD=6;
			break;;
		Change_Option )
			# if no gedit, use nano
			gedit $path2config/config.txt
			exit;;
		EXIT )	echo 'Now exiting'; 
			exit;;
	esac
done

$path2run/run.sh $SD

