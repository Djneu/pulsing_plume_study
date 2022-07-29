Files used in the study:
 
Temporal variations in plume flux: Characterizing pulsations from tilted plume conduits in a rheologically complex mantle
Derek Neuharth and Eric Mittelstaedt

Models were run using: 
deal.II 8.5.1 
ASPECT 2.0.0-pre commit 179d6da

the specific branch to run these files are at:
https://github.com/Djneu/aspect/tree/pt_test


The files contained within here consist are in 3 separate folders following the name scheme:

"rA_tzBC_vD_tE_pF"

where

A = rheology type. 1 for linear, and 2 for composite.

B = 410-km clapeyron slope.

C = 660-km clapeyron slope.

D = plate velocity

E = initial plume temperature anomaly

F = Factor that dislocation creep prefactor was multiplied by (e.g., 10 for 10x increase, and 01 for 10x decrease).

1. prmfiles
	- This contains prm files to rerun the models

2. scripts
	- This contains all the matlab scripts used to calculate everything mentioned in the text.
		  The paraview .pvsm file used to extract data is included.
          After loading this, the point and field data need to be extracted from the model using paraview.

3. numerics
	- This folder contains the statistic files, as well as all the .mat files created from the 
	  scripts in the script folder.
