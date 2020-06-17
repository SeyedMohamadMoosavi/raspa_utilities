#! /usr/bin/env python

import sys

with open ("fnames.txt") as fi:
    names = [n.strip().split()[0] for n in fi.readlines()]

KH=[]
HOA=[]
with open("kh_CO2.csv","w") as fo:
    fo.write("MOFname,filename,kH,HOA\n")
    for name in names:
        with open(name) as fi:
            data = fi.readlines()
            k=0
            h=0
            cond1 = False
            cond2 = False
            for line in data:
                if " Average Henry coefficient:" in line:
                    k = line.strip().split()[4]
                    cond1 = True
                elif " Average  <U_gh>_1-<U_h>_0:" in line:
                    h = line.strip().split()[8]
                    cond2 = True
            mofname = "_".join(st for st in name.split("_")[1:-3])
            if cond1 and cond2:
                fo.write("%s,%s,%s,%s\n"%(mofname,name,k,h))
                KH.append(k)
                HOA.append(h)
            else:
                print(name)
    
    
    
