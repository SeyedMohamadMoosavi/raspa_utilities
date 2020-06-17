#! /usr/bin/env python

## script for high-throuput analysing RASPA output files ##
##
##

import json
import numpy as np 
import sys
import math
import os
import subprocess
import matplotlib.pyplot as plt

# usage ./PP_full_isotherm workdirectory outname
# e.g. ./PP_full_isotherm.py . isotherm.csv


p=subprocess.Popen(["find",sys.argv[1],"-name","output_*"],
       stderr=subprocess.PIPE,stdout=subprocess.PIPE)
(directories , stderr)= p.communicate()
raspa_files=directories.decode("utf-8").split("\n")
raspa_files = [n for n in raspa_files if not n==""]

MOFnames = ["_".join(s for s in n.strip().split("_")[1:-3]) for n in raspa_files]
MOFnames = [n.split("output_")[-1] for n in MOFnames]
pressures = [float(n.strip().split("_")[-1].split(".")[0]) for n in raspa_files]
temperatures = [float(n.strip().split("_")[-2]) for n in raspa_files]
print(len(MOFnames))
all_data = {}
# 0. make new/update key for a MOFname
# 1. Read raspa file and find out how many components do you have
# 2. Read the adsorption data for each component and store them with proper key in the python dict
# 3. store all as a json/csv file

for raspafile,name,sim_temp,sim_p in zip(raspa_files,MOFnames,temperatures,pressures):
    with open(raspafile) as fi:
        isotherm_data = fi.readlines()
        if not isotherm_data[-2].startswith("The end time was "):
            continue
        
        file_sections = {}
        for k1,l1 in enumerate(isotherm_data):
            if set(l1.strip()) == set("="):
                for k2,l2 in enumerate(isotherm_data[k1+1:]):
                    if set(l2.strip()) == set("="):
                        break

                file_sections["%s"%isotherm_data[k1-1].strip().replace(":","")]=[k1,k1+k2]

        adsorption_data = {}
        # finding the components of the adsorption:
        k1,k2 = file_sections["MoleculeDefinitions"]
        adsorbates = []
        molfracs = []
        for k,line in enumerate(isotherm_data[k1:k2]):
            if "Adsorbate molecule" in line:
                adsorbates.append( line.strip().split("[")[1].split("]")[0])
            if "MolFraction" in line:
                molfracs.append(float(line.strip().split()[1]))

        if len(molfracs) == len(adsorbates):
            for gas,mf in zip(adsorbates,molfracs):
                adsorption_data[gas]={"molfraction":mf}
        else:
            print("Num. components doesn't match with molfracs, check this file: %s"%raspafile)
            continue

        # parsing the adsorption loadings or henry coefficients:
        if sim_p == 0:
            kHs=[]
            HOAs=[]
            k1,k2 = file_sections["Average Henry coefficient"]
            for gas in adsorbates:
                for k,line in enumerate(isotherm_data[k1:k2]):
                    if "[%s] Average Henry coefficient"%(gas) in line:
                        kHs.append( line.strip().split()[4])

            k1,k2 = file_sections["(Note the total heat of adsorption is dH=<U_gh>_1-<U_h>_0 - <U_g> - RT)"]
            for gas in adsorbates:
                for k,line in enumerate(isotherm_data[k1:k2]):
                    if " Average  <U_gh>_1-<U_h>_0:" in line:
                        HOAs.append(line.strip().split()[8])
    
            if len(kHs) == len(adsorbates):
                for gas,kH in zip(adsorbates,kHs):
                    adsorption_data[gas].update({"%s_kH"%(gas):kH})
            else:
                print("Num. components doesn't match with kH, check this file: %s"%raspafile)
                continue
    
            if len(HOAs) == len(adsorbates):
                for gas,hoa in zip(adsorbates,HOAs):
                    adsorption_data[gas].update({"%s_widomHOA"%(gas):hoa})
            else:
                print(HOAs)
                print(adsorbates)
                print("Num. components doesn't match with HOA, check this file: %s"%raspafile)
                continue
                


        else:
            k1,k2 = file_sections["Number of molecules"]
            loadings = []
            for k,line in enumerate(isotherm_data[k1:k2]):
                if "Average loading absolute [mol/kg framework]" in line:
                    loadings.append(float(line.strip().split()[5]))

            if len(loadings) == len(adsorbates):
                for gas,loading,mf in zip(adsorbates,loadings,molfracs):
                    if mf == 1:
                        adsorption_data[gas].update({"uptake_%s_pure_%.2f_%.0f"%(gas,sim_temp,sim_p):loading})
                    else:                                       
                        adsorption_data[gas].update({"uptake_%s_mix_%.2f_%.0f"%(gas,sim_temp,sim_p):loading})
            else:
                print("Num. components doesn't match with loadings, check this file: %s"%raspafile)
                continue
        
        print(adsorption_data)

    # do update or add new entry
    try:
        if len(adsorbates)>1:
            mix_name = "mix_"+"_".join('%s_%0.4f'%(gas,mf) for gas,mf in zip(adsorbates,molfracs))+"_"
            all_data[name][mix_name]=adsorption_data
        else:
            MOFname = all_data[name]
            try:
                all_data[name]["pure"][adsorbates[0]].update(adsorption_data[adsorbates[0]])
            except KeyError:
                all_data[name]["pure"]=adsorption_data
    except KeyError:
        all_data[name]={}
        if len(adsorbates)>1:
            mix_name = "mix_"+"_".join('%s_%0.4f'%(gas,mf) for gas,mf in zip(adsorbates,molfracs))+"_"
            all_data[name][mix_name]=adsorption_data
        else:
            all_data[name]["pure"]=adsorption_data
        
    
print(all_data)    
                   
with open(sys.argv[2],"w") as f:
    json.dump(all_data,f,indent=4)                   
    
