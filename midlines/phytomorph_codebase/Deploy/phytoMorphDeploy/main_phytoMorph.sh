#!/bin/bash
echo "Which program do you want to configure or EXIT?"
select selection in "phytoG" "yart" "EXIT"; do
	case $selection in
		phytoG )
			progName="phytoG"; 
			REPO="https://bitbucket.org/leeatuw/repository_experiment/raw/Jun202016";
			MCRver="v717";
			MRCticN="HlB1NJlWX8f8mqJ";
			mainP="/opt/imagePhenomics";
			MCRdownP=$mainP"/common"
			MCRruntimeVer="R2012a"
			MCRuzipP="$MCRdownP/MATLAB/MATLAB_Compiler_Runtime_$MCRruntimeVer";
			break;;
		yart ) 
			progName="mainTrace_GUI";
			REPO="https://bitbucket.org/leeatuw/repository_experiment/raw/Jun202016";
			MCRver="v717";
			MRCticN="HlB1NJlWX8f8mqJ";
			mainP="/opt/imagePhenomics";
			MCRdownP=$mainP"/common"
			MCRruntimeVer="R2012a"
			MCRuzipP="$MCRdownP/MATLAB/MATLAB_Compiler_Runtime_$MCRruntimeVer";
			break;;

		EXIT ) echo 'Now exiting'; exit;;
	esac
done

echo "Start download"
wget -O sub_phytoMorph.sh $REPO/sub_phytoMorph.sh
echo "End download"
echo "Start run "
sudo chmod +x sub_phytoMorph.sh
sudo ./sub_phytoMorph.sh $progName $REPO $MCRver $MRCticN $mainP $MCRdownP $MCRruntimeVer $MCRuzipP
echo "End run "



