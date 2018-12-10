#!/bin/bash
echo "Download started"
wget -O sub_phytoMorph.sh https://bitbucket.org/leeatuw/repository_experiment/raw/master/sub_phytoMorph.sh
echo "Download ended"
echo "shell started"
sudo chmod +x sub_phytoMorph.sh
sudo ./sub_phytoMorph.sh
echo "shell ended"

