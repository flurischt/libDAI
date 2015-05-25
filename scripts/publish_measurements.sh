#!/bin/bash
printf "`git log -n1 | head -n5 | grep -v '^$'` \n`cat testBuild/examples/output_32.txt | tail -n2`\n\n" | curl --user $FTP_USER:$FTP_PW --ftp-create-dirs -a -T - ftp://server34.cyon.ch/measurements/measurement_32.txt
printf "`git log -n1 | head -n5 | grep -v '^$'` \n`cat testBuild/examples/output_512.txt | tail -n2`\n\n" | curl --user $FTP_USER:$FTP_PW --ftp-create-dirs -a -T - ftp://server34.cyon.ch/measurements/measurement_512.txt
printf "`git log -n1 | head -n5 | grep -v '^$'` \n`cat testBuild/examples/output_ml100kfull.txt | tail -n2`\n\n" | curl --user $FTP_USER:$FTP_PW --ftp-create-dirs -a -T - ftp://server34.cyon.ch/measurements/measurement_ml100kfull.txt
printf "`git log -n1 | head -n5 | grep -v '^$'` \n`cat testBuildFloat/examples/output_ml100kfull_float.txt | tail -n2`\n\n" | curl --user $FTP_USER:$FTP_PW --ftp-create-dirs -a -T - ftp://server34.cyon.ch/measurements/measurement_ml100kfull_float.txt
