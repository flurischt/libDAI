#!/bin/bash

# $1 = cycles
# $2 = seconds
printf "`git log | head -n3` \n$1 cycles \t $2 seconds\n\n" | curl --user $FTP_USER:$FTP_PW --ftp-create-dirs -a -T - ftp://server34.cyon.ch/are_we_fast_yet.txt
