#!/bin/bash
printf "`git log -n1 | head -n5 | grep -v '^$'` \n`cat testBuild/examples/output.txt | tail -n4`\n\n" | curl --user $FTP_USER:$FTP_PW --ftp-create-dirs -a -T - ftp://server34.cyon.ch/are_we_fast_yet.txt
