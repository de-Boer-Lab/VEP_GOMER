#!/bin/bash
tar -xf CIS-BP_1.02.tar.gz
cat CIS-BP_1.02/Homo_sapiens/VEP.CISBP.BestCD.20180802.cutoffs.txt | sed 's#INSTALL_PATH#'`pwd`'#g' > Homo_sapiens.VEP.CISBP.txt
cat CIS-BP_1.02/Homo_sapiens/TF_Information.BestCD.txt  | sed 's#INSTALL_PATH#'`pwd`'#g' > Homo_sapiens.TF_Information.txt 
