#!/bin/bash
printf "`git log | head -n3` \n`cat testBuild/examples/output.txt | tail -n2`\n\n" | curl --user $FTP_USER:$FTP_PW --ftp-create-dirs -a -T - ftp://server34.cyon.ch/are_we_fast_yet.txt
