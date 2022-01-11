Files used in the study:
 
Temporal variations in plume flux: Characterizing pulsations from tilted plume conduits in a rheologically complex mantle
Derek Neuharth and Eric Mittelstaedt

Models were run using: 
deal.II 8.5.1 
ASPECT 2.0.0-pre commit 179d6da

the specific branch to run these files are at:
https://github.com/Djneu/aspect/tree/pt_test


The files contained within here consist are in 3 separate folders

1. prmfiles
	- This contains prm files to rerun the models divided into diffusion and composite
	  rheology.
	  NOTE: For the diffusion runs you may need to remove the phase tracker visualization.

2. scripts
	- This contains all the matlab scripts used to calculate everything mentioned in the text.
          To run these point and field data will need to be extracted from the model using paraview.

3. numerics
	- This folder contains the statistic files, as well as all the .mat files created from the 
	  scripts in the script folder.