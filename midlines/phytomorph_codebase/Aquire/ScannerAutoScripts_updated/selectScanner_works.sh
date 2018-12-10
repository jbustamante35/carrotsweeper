#!/bin/bash
echo "Which scanner do you want to run or EXIT?"
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
			gedit /home/spaldinglab/Desktop/production_scripts/Scripts_main/Sam_config.txt
			exit;;
		EXIT )	echo 'Now exiting'; 
			exit;;
	esac
done

/home/spaldinglab/Sam_run.sh $SD
