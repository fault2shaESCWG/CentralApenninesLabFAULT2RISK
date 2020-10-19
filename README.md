# CentralApennines
This repository contains all the files and codes necessary to run the Central Apennines model published in Scotti et al. (2020)

README
This README will guide users through the different computational steps used in Scotti et al. (2020). Please note that SHERIFS input file (FAULT2RISK /A_SHERIFS_CAD/input/ CAD_optionA1B1C1/run.info) is here parametrized to run only 1 montecarlo exploration of the multi-fault rupture model. 
We strongly suggest to perform a first calculation with the given parameters and then, if needed, to modify montecarlo exploration in the run.info file.
Data used in Scotti et al. (2020) are from the FAULT2SHA Central Apennines Database of Faure Walker et al. (2020), however, they have been also included this repository. Please note that here we only use the OPTION A1B1C1. The interested reader can re-run with alternative options, if necessary.

----------------------------------------------------------------------------------------------- 
Citation 
Complete or partial use of all the matlab scripts is allowed and has to be cited as: Scotti et al. (2020).
For SHERIFS and OPENQUAKE please cite their respective original manuscripts by Chartier et sl. (2019) and Pagani et al. (2014).
Complete or partial use of the fault data is allowed and has to be cited as: Faure Walker et al. (2020).

-----------------------------------------------------------------------------------------------
Software requirements:  
1.	SHERIFS V1.2 (Python 3.7) - available in the supplementary material with the mfd_shape.py subroutine adapted for the CAD (Mrupt = 6.7). We suggest to run test_SHERIFS.py before running SHERIFS. Please refer to the User_Manual_SHERIFS_V1.1.pdf contained in the folder A_SHERIFS_CAD for details of the code, if necessary.
2.	OPENQUAKE 3.9 (Python 3.6) - available at https://github.com/gem/oq-engine/#openquake-engine. Input and output format of the openquake files can change in different versions, we refer to the OQ 3.9 version for the file formats.
3.	MATLAB (R2016a - R2018a with mapping toolbox) - licence required
-----------------------------------------------------------------------------------------------

Steps of the calculation:
1° Running SHERIFS
2° Running OPENQUAKE
3° Visualizing

-----------------------------------------------------------------------------------------------
1° Running SHERIFS:
a.	Open Matlab and move to the folder FAULT2RISK
b.	Run A1_script_DB2SHERIFSinputs.m that will read and extract information from the DATABASE and prepare slip rate profiles, sections parameters and input files for SHERIFS (it can take several minutes). At this step you have built a fault model based on OPTIONS A1B1C1. If you want to change options you need to go back to the Excel File of the database.
c.	Run A2_script_Combine_Sections.m that will prepare a rupture list input file for SHERIFS (it can take several minutes). At this stage you have a fault model and possible rupture models based on given distance criteria and sections lengths.
d.	Open a Terminal, move to the folder “FAULT2RISK/A_SHERIFS_CAD” and type command line : 
e.	python 1_SHERIFS.py   (this program can take more than 15 minutes).
f.	python 2_Visualisation.py 
At this stage you have explored magnitude frequency distributions for each sections of your fault system and for each rupture and you have created a earthquake rupture forecast model.

2° Running hazard and risk:
g.	Open Matlab and move to the folder FAULT2RISK
h.	Run B1_script_BuildSourceModelForEachScenario.m that will prepare input files for Openquake
i.	Open a Terminal and type : source ~/openquake/env.sh
j.	From the Terminal, move to the folder FAULT2RISK and type : 
oq engine --run B_OQ_CALCULATION_GMPE_FRAGILTY_EXPOSURE/job_damage.ini --log-level info
(this run can take more than 30 minutes)
To export hazard and risk results:
k.	From the terminal type : mkdir WORKING_DIRECTORY_A1B1C1_10km/OQoutputs
l.	From the terminal type : 
oq engine --export-outputs CALCULATION/WORKING_DIRECTORY_A1B1C1_10km/OQoutputs/
NB: calculationNumber depends on your own computer. You can see this number during the OpenQuake execution. For example “calculation #2 completed in 2178 seconds”, calculationNumber = 2

3° Visualizing results
m.	Open Matlab and move to the folder FAULT2RISK
n.	Set the variable “OQ_RUN_ID” with the calculationNumber in the “user options” section of matlab scripts 

  •	C1_hazard_maps.m, 
  
  •	C2_risk_map, 
  
  •	D1_PartecipationHazard_by_section 
  
  •	D2_PartecipationRisk_by_section  
  
o.	Run the Matlab codes
p.	Figures are in WORKING_DIRECTORY_A1B1C1_10km/visualization/figures
