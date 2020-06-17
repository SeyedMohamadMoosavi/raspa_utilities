#! /usr/bin/env python

import numpy as np
import sys
import subprocess
import math
from atomic import UFF_atoms



uff_atoms = [a.lower() for a  in UFF_atoms]

### functions ###
deg2rad = math.pi/180.0
########## Functions ################# 
def writecif(file_name,cell,atoms):
    with open(file_name , 'w') as fo:
        ## writing the general P1 symmetry ##
        fo.write("data_%s\n"%file_name)
        fo.write("_audit_creation_method           standardcif\n")
        fo.write("_symmetry_space_group_name_H-M    P1\n")
        fo.write("_symmetry_Int_Tables_number       1\n")
        fo.write("_symmetry_cell_setting            triclinic\n")
        fo.write("loop_\n")
        fo.write("_symmetry_equiv_pos_as_xyz\n")
        fo.write("\'x, y, z\'\n")
        cell_a=math.sqrt(cell[0][0]**2+cell[0][1]**2+cell[0][2]**2)
        cell_b=math.sqrt(cell[1][0]**2+cell[1][1]**2+cell[1][2]**2)
        cell_c=math.sqrt(cell[2][0]**2+cell[2][1]**2+cell[2][2]**2)
        cell_gamma=math.acos(cell[1][0]/cell_b)
        cell_beta =math.acos(cell[2][0]/cell_c)
        cell_alpha=math.acos(cell[2][1]/cell_c*math.sin(cell_gamma)+math.cos(cell_beta)*math.cos(cell_gamma))
        cell_gamma=np.degrees( cell_gamma)
        cell_beta =np.degrees( cell_beta )
        cell_alpha=np.degrees( cell_alpha)
        fo.write("_cell_length_a %12.08f\n"%(cell_a))
        fo.write("_cell_length_b %12.08f\n"%(cell_b))
        fo.write("_cell_length_c %12.08f\n"%(cell_c))
        fo.write("_cell_angle_alpha %12.08f\n"%(cell_alpha))
        fo.write("_cell_angle_beta  %12.08f\n"%(cell_beta ))
        fo.write("_cell_angle_gamma %12.08f\n"%(cell_gamma))
        fo.write("loop_\n")
        fo.write("_atom_site_type_symbol\n")
        fo.write("_atom_site_fract_x\n")
        fo.write("_atom_site_fract_y\n")
        fo.write("_atom_site_fract_z\n")
        fo.write("_atom_site_label\n")
        charge_flag = False
        if "charge" in atoms[0].keys():
            fo.write("_atom_site_charge\n")
            charge_flag = True
        for cn,atom in enumerate(atoms):
            if charge_flag:
                fo.write("%-10s %10.6f %10.6f %10.6f %10s %10.6f\n"%(
                    atom["atom_type"]
                    ,atom["fracx"]
                    ,atom["fracy"]
                    ,atom["fracz"]
                    ,atom["atom_type"]
                    #,atom["atom_id"]
                    ,atom["charge"]))
            else:
                fo.write("%-10s %10.6f %10.6f %10.6f %10s\n"%(
                    atom["atom_type"]
                    ,atom["fracx"]
                    ,atom["fracy"]
                    ,atom["fracz"]
                    ,atom["atom_type"]
                    ))
    
def readcif(name):
    with open (name , 'r') as fi:
        EIF = fi.readlines()
        cond=True
        cond2=False
        atom_props_count=0
        atomlines=[]
        counter=0
        cell_parameter_boundary=[0.0,0.0]
        for line in EIF:
            line_stripped=line.strip()
            if (not line) or line_stripped.startswith("#"):
                continue
            line_splitted=line.split()
            if line_stripped.startswith("_cell_length_a"):
                cell_a=line_splitted[1]
                cell_parameter_boundary[0]=counter+1
            elif line_stripped.startswith("_cell_length_b"):
                cell_b=line_splitted[1]
            elif line_stripped.startswith("_cell_length_c"):
                cell_c=line_splitted[1]
            elif line_stripped.startswith("_cell_angle_alpha"):
                cell_alpha=line_splitted[1]
            elif line_stripped.startswith("_cell_angle_beta"):
                cell_beta=line_splitted[1]
            elif line_stripped.startswith("_cell_angle_gamma"):
                cell_gamma=line_splitted[1]
                cell_parameter_boundary[1]=counter+1
            if cond2==True and line_stripped.startswith("loop_"):
                break
            else:
                if line_stripped.startswith("_atom") :
                    atom_props_count+=1
                    if line_stripped=="_atom_site_label":
                        type_index=atom_props_count-1
                    elif line_stripped=="_atom_site_fract_x":
                        fracx_index=atom_props_count-1
                    elif line_stripped=="_atom_site_fract_y":
                        fracy_index=atom_props_count-1
                    elif line_stripped=="_atom_site_fract_z":
                        fracz_index=atom_props_count-1
                    elif "charge" in line_stripped:
                        charge_index=atom_props_count-1
        
                    cond2=True
                elif cond2==True:
                    if len(line_splitted)==atom_props_count:
                        atomlines.append(line)
        
            counter+=1
        
        positions=[]
        numbers=[]
        types=[]
        atoms=[]
        for count,at in enumerate(atomlines):
            ln=at.strip().split()
            positions.append([float(ln[fracx_index]),float(ln[fracy_index]),float(ln[fracz_index])])
            types.append(ln[type_index])
            try:
                atoms.append({"atom_id":count,"atom_type":''.join([iii for iii in ln[type_index] if not iii.isdigit()]),
                        "charge":float(ln[charge_index])
                       ,"fracx" :float(ln[fracx_index])
                       ,"fracy" :float(ln[fracy_index])
                       ,"fracz" :float(ln[fracz_index])})
            except NameError:
                atoms.append({"atom_id":count,"atom_type":''.join([iii for iii in ln[type_index] if not iii.isdigit()])
                                ,"fracx":float(ln[fracx_index])
                                ,"fracy":float(ln[fracy_index])
                                ,"fracz":float(ln[fracz_index])})

            
        unique_types=list(set(types))
        
        for label in types:
            numbers.append(int(unique_types.index(label)+1))
        cpar=[cell_a,cell_b,cell_c]
        cang=[cell_alpha,cell_beta,cell_gamma]
        cpar=[float(i) for i in cpar]
        cang=[float(i)*deg2rad for i in cang]
        v=math.sqrt(1-math.cos(cang[0])**2-math.cos(cang[1])**2-math.cos(cang[2])**2+2*math.cos(cang[0])*math.cos(cang[1])*math.cos(cang[2]))
        CP=[]
        CP.append([cpar[0],0,0])
        CP.append([cpar[1]*math.cos(cang[2]),cpar[1]*math.sin(cang[2]),0])
        CP.append([cpar[2]*math.cos(cang[1]),cpar[2]*(math.cos(cang[0])-math.cos(cang[1])*math.cos(cang[2]))/(math.sin(cang[2])),cpar[2]*v/math.sin(cang[2])])
        CP=np.array(CP)
        return CP, atoms


# def find_min_cell(cell):
#     ew;lgj


### Simulation variables ###




class RASPA_calculation(object):
    def __init__(self):
        self.IniCycles= 2000
        self.TotCycles= 6000
        self.SimulationType= "MonteCarlo"
        self.PrintEvery = 1000
        self.PrintPropertiesEvery=1000
        self.ContinueAfterCrash="no"
        self.RestartFile="no"
        self.RestartStyle="Raspa"
        self.CutOff = 12.8
        self.Forcefield = "Local"
        self.FrameworkName = None 
        self.InputFileType = "cif"
        self.UseChargesFromCIFFile = "yes"   
        self.UnitCells = [1,1,1]              
        self.HeliumVoidFraction = 0.0     
        self.Movies = "yes"                  
        self.WriteMoviesEvery = 100       
        self.ExternalTemperatures=None
        self.ExternalPressures=None
        self.gascomposition = None
        self.cifnames = None 
        self.path2cifs=None
        self.calcname=None
        self.spacingVDWG=0.5
        self.spacingCoulombG=0.5
        self.MoleculeDefinition="local"

    def setup(self,names,p2c,calctype,calcname):
        self.cifnames = names
        self.path2cifs = p2c
        self.calcname=calcname
        if calctype =="GCMC":
            self.SimulationType= "MonteCarlo"
        elif calctype == "egrid":
            self.SimulationType = "MakeEGrid"
        elif calctype == "henry":
            self.SimulationType = "henrycoeff"
        else:
            print("%s type of calculation is not supported..\n "%calctype)
            sys.exit()
        
    def compute_min_rep(self,cell):
        """Calculate the smallest supercell with a half-cell width cutoff.

        Increment from smallest cell vector to largest. So the supercell
        is not considering the 'unit cell' for each cell dimension.

        """

        a_cross_b = np.cross(cell[0], cell[1])
        b_cross_c = np.cross(cell[1], cell[2])
        c_cross_a = np.cross(cell[2], cell[0])

        #volume = np.dot(self.cell[0], b_cross_c)

        widths = [np.dot(cell[0], b_cross_c) / np.linalg.norm(b_cross_c),
                  np.dot(cell[1], c_cross_a) / np.linalg.norm(c_cross_a),
                  np.dot(cell[2], a_cross_b) / np.linalg.norm(a_cross_b)]

        return tuple(int(math.ceil(2*self.CutOff/x)) for x in widths)

    def write_sim_input(self,n,p,t):
        c,ats = readcif(self.path2cifs+self.cifnames[n]+".cif")
        for at in ats:
            if not at["atom_type"].lower() in uff_atoms:

                print("\nStructure %s has atom type %s which is not compatible with UFF\n"%(self.cifnames[n],at["atom_type"].lower()))
                with open("failed_structures.txt","a") as f_fail:
                    f_fail.write("%s,structure %s has atom type %s which is not compatible with UFF\n"%(self.cifnames[n],self.cifnames[n],at["atom_type"].lower()))

                return 0 


        writecif(self.cifnames[n]+".cif",c,ats)

        if self.SimulationType== "MonteCarlo":
            with open("%s.input"%self.cifnames[n],"w") as fo:
                ## simulation part ##
                str_out = ""
                str_out += "SimulationType               %s\n"%(self.SimulationType) 
                str_out += "NumberOfCycles               %i\n"%(self.TotCycles) 
                str_out += "NumberOfInitializationCycles %i\n"%(self.IniCycles) 
                str_out += "PrintEvery                   %i\n"%(self.PrintEvery) 
                str_out += "PrintPropertiesEvery         %i\n"%(self.PrintPropertiesEvery) 
                str_out += "ContinueAfterCrash           %s\n"%(self.ContinueAfterCrash) 
                str_out += "RestartFile                  %s\n"%(self.RestartFile) 
                str_out += "RestartStyle                 %s\n"%(self.RestartStyle) 
                str_out += "CutOff                       %4.3f\n"%(self.CutOff) 
                str_out += "Forcefield                   %s\n"%(self.Forcefield) 
                str_out += "Framework 0                  \n\n"
                str_out += "FrameworkName                %s\n"%(self.cifnames[n]) 
                str_out += "InputFileType                %s\n"%(self.InputFileType) 
                str_out += "UseChargesFromCIFFile        %s\n"%(self.UseChargesFromCIFFile) 
                (rep1,rep2,rep3) = self.compute_min_rep(c)
                str_out += "UnitCells              %i %i %i\n"%(rep1,rep2,rep3) 
                str_out += "HeliumVoidFraction        %4.2f\n"%(self.HeliumVoidFraction) 
                str_out += "Movies                       %s\n"%(self.Movies) 
                str_out += "WriteMoviesEvery             %i\n"%(self.WriteMoviesEvery) 
                str_out += "ExternalTemperature       %8.3f\n"%(self.ExternalTemperatures[t]) 
                print(self.ExternalPressures,p)
                str_out += "ExternalPressure    %8.3f\n\n\n\n"%(self.ExternalPressures[p])
                if len(self.gascomposition.keys())==1:
                    g = list(self.gascomposition.keys())[0]
                    str_out += "Component 0 MoleculeName    %s\n"%(g)
                    str_out += "MoleculeDefinition           %s\n"%self.MoleculeDefinition 
                    str_out += "TranslationProbability      0.5\n"
#                    str_out += "RegrowProbability           0.5\n" # the keyword in raspa is different for a reason!?
                    str_out += "ReinsertionProbability      0.5\n"
                    str_out += "SwapProbability             1.0\n"
                    str_out += "CreateNumberOfMolecules       0\n"
                    str_out +="\n\n"
                else:
                    for ii,g in enumerate(self.gascomposition.keys()):
                        str_out += "Component %i MoleculeName    %s\n"%(ii,g)
                        str_out += "MoleculeDefinition           %s\n"%self.MoleculeDefinition 
                        str_out += "MolFraction               %10.08f\n"%(self.gascomposition[g])
                        str_out += "TranslationProbability      0.5\n"
                        str_out += "RegrowProbability           0.5\n"
                        str_out += "IdentityChangeProbability   1.0\n"
                        str_out += "  NumberOfIdentityChanges    %i\n"%(len(self.gascomposition.keys()))
                        tmp_str = " ".join([str(iii) for iii in range(len(self.gascomposition.keys()))])
                        str_out += "  IdentityChangesList       %s\n"%tmp_str
                        str_out += "SwapProbability             1.0\n"
                        str_out += "CreateNumberOfMolecules       0\n"
                        str_out +="\n\n"


                fo.write(str_out)

        if self.SimulationType== "henrycoeff":
            with open("%s.input"%self.cifnames[n],"w") as fo:
                ## simulation part ##
                str_out = ""
                str_out += "SimulationType               %s\n"%("MonteCarlo") 
                str_out += "NumberOfCycles               %i\n"%(self.TotCycles) 
                str_out += "NumberOfInitializationCycles %i\n"%(0) 
                str_out += "PrintEvery                   %i\n"%(self.PrintEvery) 
                str_out += "PrintPropertiesEvery         %i\n"%(self.PrintPropertiesEvery) 
                str_out += "CutOff                       %4.3f\n"%(self.CutOff) 
                str_out += "Forcefield                   %s\n"%(self.Forcefield) 
                str_out += "Framework 0                  \n\n"
                str_out += "FrameworkName                %s\n"%(self.cifnames[n]) 
                str_out += "InputFileType                %s\n"%(self.InputFileType) 
                str_out += "UseChargesFromCIFFile        %s\n"%(self.UseChargesFromCIFFile) 
                (rep1,rep2,rep3) = self.compute_min_rep(c)
                str_out += "UnitCells              %i %i %i\n"%(rep1,rep2,rep3) 
                str_out += "ExternalTemperature       %8.3f\n\n\n"%(self.ExternalTemperatures[t]) 
                if len(self.gascomposition.keys())==1:
                    g = list(self.gascomposition.keys())[0]
                    str_out += "Component 0 MoleculeName    %s\n"%(g)
                    str_out += "MoleculeDefinition           %s\n"%self.MoleculeDefinition 
                    str_out += "WidomProbability            1.0\n"
                    str_out += "CreateNumberOfMolecules       0\n"
                    str_out +="\n\n"
                else:
                    print("multicomponent henry coefficient doesn't make sense!\n")
                    sys.exit()

                fo.write(str_out)

        elif self.SimulationType== "MakeEGrid":
            with open("%s.input"%self.cifnames[n],"w") as fo:
                ## simulation part ##
                str_out = ""
                str_out += "SimulationType               %s\n\n"%(self.SimulationType) 
                str_out += "CutOff                       %4.3f\n"%(self.CutOff) 
                str_out += "Forcefield                   %s\n"%(self.Forcefield) 
                str_out += "Framework 0                  \n\n"
                str_out += "FrameworkName                %s\n"%(self.cifnames[n]) 
                str_out += "InputFileType                %s\n"%(self.InputFileType) 
                str_out += "UseChargesFromCIFFile        %s\n"%(self.UseChargesFromCIFFile) 
                (rep1,rep2,rep3) = self.compute_min_rep(c)
                str_out += "UnitCells              %i %i %i\n"%(rep1,rep2,rep3) 
                if len(self.gascomposition.keys())==1:
                    g = list(self.gascomposition.keys())[0]
                    str_out += "NumberOfGrids  %i\n"%(len(self.gascomposition.keys()))
                    str_out += "GridTypes    %s\n"%(g)
                    str_out += "SpacingVDWGrid       %f\n"%self.spacingVDWG
                # else:
                #     for ii,g in enumerate(self.gascomposition.keys()):
                #         str_out += "Component %i MoleculeName    %s\n"%(ii,g)
                #         str_out += "MoleculeDefinition           %s\n"%self.MoleculeDefinition 
                #         str_out += "MolFraction               %10.08f\n"%(self.gascomposition[g])
                #         str_out += "TranslationProbability      0.5\n"
                #         str_out += "RegrowProbability           0.5\n"
                #         str_out += "IdentityChangeProbability   1.0\n"
                #         str_out += "  NumberOfIdentityChanges    %i\n"%(len(self.gascomposition.keys()))
                #         tmp_str = " ".join([str(iii) for iii in range(len(self.gascomposition.keys()))])
                #         str_out += "  IdentityChangesList       %s\n"%tmp_str
                #         str_out += "SwapProbability             1.0\n"
                #         str_out += "CreateNumberOfMolecules       0\n"
                #         str_out +="\n\n"
                fo.write(str_out)
        return 1

    def write_slurm_job(self,n):
        if self.SimulationType in ["MonteCarlo","henrycoeff"]:
            str_out = ""
            str_out += "#!/bin/bash\n"
            str_out += "##SBATCH --partition debug\n" 
            str_out += "#SBATCH --time 72:00:00\n"
            str_out += "#SBATCH --nodes  1\n" 
            str_out += "#SBATCH --ntasks 1\n\n"

            str_out += "srun ~/RASPA_bin/bin/simulate simulation.input"
            with open ("job_%s"%n,"w") as fj:
                fj.write(str_out)
        elif self.SimulationType== "MakeEGrid":
            str_out = ""
            str_out += "#!/bin/bash\n"
            str_out += "##SBATCH --partition debug\n" 
            str_out += "#SBATCH --time 72:00:00\n"
            str_out += "#SBATCH --workdir ./\n"
            str_out += "#SBATCH --nodes  1\n" 
            str_out += "#SBATCH --ntasks 1\n\n"

            str_out += "srun ~/RASPA_mod/bin/simulate simulation.input"
            with open ("job_%s"%n,"w") as fj:
                fj.write(str_out)

    def prepare_simulations(self):
        p = subprocess.Popen(["mkdir",self.calcname])
        if self.SimulationType in ["MonteCarlo"]:
            for pid,pressure in enumerate(self.ExternalPressures):
                for tid,temperature in enumerate(self.ExternalTemperatures):
                    for nid,name in enumerate(self.cifnames):
                        simfiles_done = self.write_sim_input(nid,pid,tid)
                        if simfiles_done:
                            self.write_slurm_job(name)
                            strdir = "%s/%s"%(self.calcname,name)
                            # print(strdir)
                            wkdir = "%s/%s/%3.0f_%08.0f"%(self.calcname,name,temperature,pressure)
                            p = subprocess.Popen(["mkdir",strdir]).wait()
                            p = subprocess.Popen(["mkdir",wkdir]).wait()
                            p = subprocess.Popen(["mv","%s.cif"%name,wkdir]).wait()
                            p = subprocess.Popen(["mv","%s.input"%name,"%s/simulation.input"%wkdir]).wait()
                            p = subprocess.Popen(["mv","job_%s"%name,"%s/job_%s"%(wkdir,name)]).wait()
                            p = subprocess.Popen(["cp","pseudo_atoms.def","force_field_mixing_rules.def"
                                                    ,"%s/"%(wkdir)]).wait()
                            # p = subprocess.Popen(["cp","force_field.def"
                            #                         ,"%s/"%(wkdir)]).wait()
                            for g in self.gascomposition.keys():
                                p = subprocess.Popen(["cp","%s.def"%g
                                                                        ,"%s/"%(wkdir)]).wait()
            
        if self.SimulationType in ["henrycoeff"]:
            for tid,temperature in enumerate(self.ExternalTemperatures):
                for nid,name in enumerate(self.cifnames):
                    simfiles_done = self.write_sim_input(nid,0,tid)
                    if simfiles_done:
                        self.write_slurm_job(name)
                        strdir = "%s/%s"%(self.calcname,name)
                        print(strdir)
                        wkdir = "%s/%s/%3.0f"%(self.calcname,name,temperature)
                        p = subprocess.Popen(["mkdir",strdir]).wait()
                        p = subprocess.Popen(["mkdir",wkdir]).wait()
                        p = subprocess.Popen(["mv","%s.cif"%name,wkdir]).wait()
                        p = subprocess.Popen(["mv","%s.input"%name,"%s/simulation.input"%wkdir]).wait()
                        p = subprocess.Popen(["mv","job_%s"%name,"%s/job_%s"%(wkdir,name)]).wait()
                        p = subprocess.Popen(["cp","pseudo_atoms.def","force_field_mixing_rules.def"
                                                ,"%s/"%(wkdir)]).wait()
                        p = subprocess.Popen(["cp","force_field.def"
                                                ,"%s/"%(wkdir)]).wait()
                        for g in self.gascomposition.keys():
                            p = subprocess.Popen(["cp","%s.def"%g
                                                                    ,"%s/"%(wkdir)]).wait()
                
        elif self.SimulationType == "MakeEGrid":
            for nid,name in enumerate(self.cifnames):
                simfiles_done = self.write_sim_input(nid,0,0)
                if simfiles_done:
                    self.write_slurm_job(name)
                    strdir = "%s/%s"%(self.calcname,name)
                    print(strdir)
                    wkdir = "%s/%s"%(self.calcname,name)
                    p = subprocess.Popen(["mkdir",strdir]).wait()
#                    p = subprocess.Popen(["mkdir",wkdir]).wait()
                    p = subprocess.Popen(["mv","%s.cif"%name,wkdir]).wait()
                    p = subprocess.Popen(["mv","%s.input"%name,"%s/simulation.input"%wkdir]).wait()
                    p = subprocess.Popen(["mv","job_%s"%name,"%s/job_%s"%(wkdir,name)]).wait()
                    p = subprocess.Popen(["cp","pseudo_atoms.def","force_field_mixing_rules.def"
                                            ,"%s/"%(wkdir)]).wait()
#                 p = subprocess.Popen(["cp","force_field.def"
#                                         ,"%s/"%(wkdir)]).wait()
#                for g in self.gascomposition.keys():
#                    p = subprocess.Popen(["cp","%s.def"%g
#                                                            ,"%s/"%(wkdir)]).wait()

        

if __name__ == "__main__":
    # run the program with:
    # run_raspa.py filenames cif_dir flag
    # example: ../run_raspa.py filenames cifs/ egrid/GCMC/henry flag 
    # example: ./run_raspa.py filenames cifs/ egrid test 
    # TODO add a check if the files required do not exist
    with open(sys.argv[1],"r") as fnames:
        cifnames=[n.strip().replace(".cif","") for n in fnames.readlines()]

    # cell,atoms=readcif("%s%s"%("PMOF_verify/cifs/",cifnames[0]))
    # writecif("test.cif",cell,atoms)
    sim = RASPA_calculation()
    sim.setup (cifnames,sys.argv[2] , sys.argv[3] ,sys.argv[4])
    # sim.gascomposition = {"N2":0.85,"CO2":0.15}
    # sim.gascomposition = {"methane":1}
    sim.gascomposition = {"CO2":1}
    # sim.gascomposition = {"N2":1.0}
    # sim.ExternalTemperatures=[298]
    # sim.ExternalPressures=[1e5]
    # sim.ExternalPressures=[1e4]
    # sim.ExternalPressures=[1.5e4,16e5]
    # sim.ExternalPressures=[16e5]
    # sim.ExternalPressures=[1e3,2.5e3,5e3,7.5e3,1e4,2.5e4,5e4,1e5]
    sim.ExternalPressures=[1e3,2.5e3,5e3,7.5e3,1e4,2.5e4,5e4,7.5e4,1e5,1.5e5,2.5e5,5e5]
    # sim.ExternalPressures=[5.8e5,65e5]
    # sim.ExternalPressures=[5.8e5]
    # sim.gascomposition = {"N2":1}
    # sim.ExternalTemperatures=[298]
    sim.ExternalTemperatures=[265,320]
    # sim.ExternalTemperatures=[363]
    # sim.ExternalPressures=[1e4]
    sim.prepare_simulations()

