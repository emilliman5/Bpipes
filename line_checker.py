#!/usr/bin/env python
import csv
import sys
import re

f = csv.reader(sys.stdin, dialect="excel-tab")
of = csv.writer(sys.stdout, dialect="excel-tab")
last_read = None
XS=0
for line in f :
    #take care of the header
    if(line[0][0] == "@") :
        of.writerow(line)
        continue

    if(last_read == None) : 
        last_read = line
        if(re.split('\t|:',line)[20] =="XS" & re.split('\t|:',line)[20] <= re.split('\t|:',line)[23]):
            XS=1
        
    else :
        if(last_read[0] == line[0]):    
            if(XS==1 & re.split('\t|:',line)[20] =="XS" & XS==0 & re.split('\t|:',line)[20] <= re.split('\t|:',line)[23]):
                last_read = None
                XS=0 
                
            else:
                 of.writerow(last_read)
                 of.writerow(line)
                 XS=0
                 last_read=None
        else :
            last_read = line
            XS=0
            if(re.split('\t|:',line)[20] =="XS" & re.split('\t|:',line)[20] <= re.split('\t|:',line)[23]):
                XS=1