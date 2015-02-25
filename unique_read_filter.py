#!/usr/bin/env python
import csv
import sys
import re

f = csv.reader(sys.stdin, dialect="excel-tab")
of = csv.writer(sys.stdout, dialect="excel-tab")

last_read = None
XS=0

for line in f:
    #take care of the header
    if(line[0][0] == "@"):
        of.writerow(line)
        continue

    if(last_read == None): 
        last_read = line
        
        if(re.split(':',line[12])[0] == "XS" and int(re.split(':',line[12])[2]) == int(re.split(':',line[12])[2])):
            XS=1
        
    else :
        if(last_read[0] == line[0]):    
            if(XS==1 and re.split(':',line[12])[0] =="XS" and int(re.split(':',line[12])[2]) == int(re.split(':',line[12])[2])):
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
            if(re.split(':',line[13])[0] == "XS" and int(re.split(':',line[12])[2]) == int(re.split(':',line[13])[2])):
                XS=1