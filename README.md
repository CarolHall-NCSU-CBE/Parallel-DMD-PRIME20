# Parallel DMD/PRIME20 - Coarse-grained force field with discontinuous molecular dynamics simulations for peptide self-assembly modelling 
## Table of contents
* [Introduction](#introduction)
* [Requirement and Installation](#requirement-and-installation)
* [Getting Started](#getting-started)
* [Running Simulation](#running-simulation)
* [Developing Status](#developing-status)
## Introduction
PRIME20 is a coarse-grained, implicit-solvent, intermediate-resolution protein model that was developed by the Hall group at North Carolina State University. The model was designed to be used with discontinuous molecular dynamics simulations (DMD) to investigate self-assembly of short peptides from their random denatured states. We are developing the parallel version of DMD/PRIME20 to reduce the computational cost for DMD simulations. PRIME20 contains geometric and energetic parameters that describe the sidechain-sidechain interactions of all 20 natural amino acids. In PRIME20, each amino acid is represented by four beads: one for the amino group (NH), one for the alpha carbon (CαH), one for the carbonyl group (CO), and one for the side chain (R). DMD/PRIME20 simulation systems are canonical ensemble (NVT) with constant number of molecules (N), simulation box volume (V) and simulation temperature (T). Temperature is maintaned by using Anderson thermostat. Neutral pH water solvent is described implicitly within the force-field. Peptides that are built and simulated by PRIME20 are capped at both terminus. DMD/PRIME20 has been used successfully to simulate spontaneous α-helix, β-sheet, and amyloid fibril formation starting from the denatured conformations of peptides such as prion proteins fragments[1][2], tau protein fragments[3], Aβ16-22 peptides[4][5][6], and  Aβ17-36 peptides[7]. DMD/PRIME20 is also used with PepAD - a peptide design software, to discover amyloid-forming peptides [8].

## Requirement and Installation
- The package has been developed since 2001 using Fortran90
- Parallelizing is done using Message Passing Interface (MPI)
- OpenMPI compiler for Fortran `mpif90` is required. OpenMPI Fortran Compiler must be installed on your device or the module that contains the compiler must be loaded before compiling. 
- The installation is through the terminal.
- The source codes are in `/src/`. To compile, open a terminal and then navigate to the `/src/` directoy on your local device. Once in '/src/' directory, create the executed files by enter the commands below.
1. To create **initconfig** for generating initial configuration
>
 	make -f genconfig.mk

2. To create **DMDPRIME20** for DMD simulations
>
	make -f dmd.mk 
- If there is no error return, check if **initconfig** and **DMDPRIME20** are succesfully created in *src*
- Obtain the paths to these executable files to use in job submission.
>**Note:** if redownload the package or update a new version, the previous steps need to be redo. 

## Getting Started   
**/example/**: this directory contains an example of required file and subdirectories for a simulation using DMD/PRIME20.
Requirements to start a simulation including **input.txt** and **parallelscript.csh** and 5 empty directories to record simulation ouputs - **/checks/**, **/inputs/**, **/outputs/**, **/parameters/**, and **/results/**.
1. **input.txt**: Please follow the format to enter all parameters that are required for a simulation. The explanation for each parameters are also included in the file.
	- **pep1** and **pep2** are sequences of the peptides that are simulated. It must be in abbrevating alphabetical format (e.g. pep1=GVLYVGS) . The current version can run simulations for system with single or double components; each with maximum length of 30 residues. If system contains single peptide sequence, then *pep1* and *pep2* must be the same in the 'input.txt'

	- **chain1** and **chain2** are the number of peptide chains of each peptide component in the simulation box. If the peptide is long, *chain1* and *chain2* should be reduced to avoid overcrowding, overlapping and to reduce simulation time. The largest system has been simulated using DMD-PRIME20 contains 200 peptides chains.

	- **boxlength** is the length of the simulation box. DMD-PRIME20 uses cubic box with periodic boundary condition for all simulations. *boxlength* is selected based on the number of peptide chains and concentration:

$$ boxlength = (\frac{\text{Total number of peptide chains}*1000}{\text{Avogadro's number * Concentration}})^\frac{1}{3}*10^9 $$
* where *Concentration* is in *mM* and *boxlength* is in *Angstrom*

	- **T** is simulation temperature in *Kelvin*. When start simulations for a new system, it is recommended to run multiple simulations of the same system at different temperatures. Check the simulation results to select the temperature that predict high order peptide aggregation. The simulation might get stuck in local miminima if the temperature is too low, but there is no aggregation if the temperature is too low.
 	- **coll** which is the number of collisions for DMD-PRIME20 to finish a *round* and record simulation results. DMD-PRIME20 is designed to run, complete and record in many rounds to avoid large result files and to allow the simulation to restart if it is crashed midway. As DMD is discontinous molecular dynamics simulation, collsion (coll) is used instead of timestep. Collision will be converted to real time when running data analysis package (underdevelopment and will be updated soon). There is not a fix value in real time for a collision.    

	- **Annealing**: The current version allows annealing simulation with a default set of temperatures (annealing = 0) or a user-defined temperatures (annealing = 1). If using user-defined temperature, include addtional parameters below the annealing line:
 		- **starting_temp**: starting temperature for the annealing process (in *Kelvin*)
		- **ending_temp**: ending temperature for the annealing process (in *Kelvin*)
		- **temp_step**: temperature drop for each annealing cycle (in *Kelvin*)
  		- **annealing_coll**: number of collisions for each temperature in annealing process. Recommended value is from 100 		million to 250 million collisions   

>**Note:** If an error is returned and the simulation is terminated during the generating of initital configuration. Adding another parameter to the very end of **input.txt**: 
>
>	*sidechainmove* = value that is larger than 3.0	
>	
>It is recommended to increase only 0.5 at a time starting from 3.0. A very large number will make the initial configuration generation very slow`

- An example of **input.txt** that include parameters for are use-defined annealing temperature is below. If running simulaiton with default annealing temperature, set annealing = 0 and delete all parameter below that line.
>
	Peptide sequence 1

	pep1=GVLYVGS
 
	Number of peptide 1 chains in the system
 
	chain1=3
 
	Peptide sequence 2
 
	pep2=GVLYVGS
 
	Number of peptide 2 chains in the system
	
 	chains=3
	
 	Box length in Angstrom

 	boxlength=159.0D0
	
 	Simulation temperature in Kelvin
	
 	T=310.0D0
	
 	Result recording frequency in collisions
	
 	coll=1000000000
	
 	Annealing process only for new simulation: 0 is the default, 1 is the specified temperatures in Kelvin
	
 	annealing = 1
	
 	max = 1000.0
	
 	min = 375.0
	
 	temp_step = 125
	
 	annealing_coll = 100000000
2. **parallelscript.csh** is an example of the tcsh script that is used to submit a job on an HPC system. This file will need to be modified according to users' computer system. Main content of the script is the three steps of a simulation.
>
	#Generate initial configuration for new simulation:

	/path_to_initconfig/initconfig

	#Annealing:

	foreach i (`seq 1 annealing_rounds`)

		mpirun /path_to_DMDPRIME20/DMDPRIME20 < inputs/annealtemp_$i > outputs/out_annealtemp_$i

	end

	#DMD simulations

	foreach i (`seq 1 simulation_rounds`)

		mpirun /path_to_DMDPRIME20/DMDPRIME20 < inputs/simtemp > outputs/out_simtemp_$i

	end

- Generating initial configuration for new simulation: This step is to create a cubic box that contents the number of peptide chains defined by users, position and velocity of each particles. Outputs of this step are saved in **/inputs/**, **/parameters/**, and **/results/** directories. These files are required for any DMD-PRIME20 simulation and need to be available in their designated locations. If restarting or resuming a simulation, this step is skipped as long as the initial configuration files are available. Although, DMD-PRIME20 simulation is parallelized, this step is done in *serial*.
- Annealing: this step is to heat up the initial system to very high temperature and then slowly cool it down to near simulation temperature. This step is only required for simulation of a completely new system. The purpose of this step is to make sure all peptide chains are denatured and simulation starts with all random coils. There are two options for annealing:
	- Default annealing (annealing = 0 in input.txt): The annealing process will be done with a default set of temperatures. These temperatures are used in many simulations since the software was developed. If using default annealing, set **annealing_rounds** in the parallelscript.sch to **9**. This means the anneanling process runs at 9 different temperatures. The temperatures and number of collision at each temperature can be found in */inputs/* directory.
 	- User-defined annealing (annealing = 1 in input.txt): The annealing process will be done with the temperature range and number of collision that are defined by user. If using this option, the **annealing_rounds** is found as:

     $$ \text{annealing_rounds} = \frac{\text{Starting temperature - Ending temperature}}{\text{temp_step}} + 1 $$




- The **.out** file shows an example of successful initial configuration generation. If your screen-written output look like this and no error showed, the initial configuration is successulffy generated. This *.out* file must be deleted before any simulation if it exists to avoid being confused by old data.  
- 5 empty directories for data recording must be created before submitting a job. The names of these directories must be exact.
	- `/checks/`: files for checking if the initial configuration is created correctly
	- `/inputs/`: files to record residue id (identity.inp and identity2.inp), positions for each peptide sequence (chninfo-n1.data and chninfo-n2.data), reduced annealing temperatures (annealtemp_*), and reduced simulation temperature (simtemp)   
	- `/outputs/`: output files for each simulation round
	- `/parameters/`: sidechain parameters generated from the inital configuration step that are required for simulation steps
	- `/results/`:  simulation results for data analysis
		1. .bptnr: collision, bond partner of each particle
		2. .config: collision, time, particle coordinates
		3. .energy: collision, time, kinetic energy, total energy, etc.
		4. .lastvel: collision, velocities 
		5. .pdb: pdb file
		6. .rca: distance from sidechain to each particle in the backbone of a residue
>Note: These subdirectories in the **/example/** directory contains results from a short simulation for your reference. When running a new simulation, these subdirectories must be empty to avoid incorrectly data appending. When running a continuing simulation, keep all results from previous simulation in these directories.


## Running simulation
DMD simulation using PRIME20 starts with building initial configuration. The current version is effective for system of no more than 30-residue peptides. It is recommended that concentration and number of peptide chains are reduced for longer peptides to avoid overlap due to overcrowding. User should check output file for overlapping error and reduce system size (number of peptides or concentration) if error is reported. PRIME20 allows simulations of a homogenous system or a heterogeneous system of two different peptides. Concentration in DMD-PRIME20 is defined by numbers of peptide chains and the simulation box length which is found as:





### Submit a job:
Steps to submit a simulation is as follow. These steps are after the package is succesfully installed on your device and *the path to executable file is obtained*.
1. Make a directory to run simulation or copy over the /'example/' directory, rename and then delete all files within subdirectories and *.out*. If making new directory, follow the next steps. 
2. In this directory, make an 'input.txt' file following the example. You can copy over this file and change the parameters correspoding to your system.
3. In this directory, make 5 empty subdirectories at listed above if running a new simulation, or copy over these subdirectories with all data in them for a continuing simulation. 
4. Submit job. It is not recommended to run DMD/PRIME20 on a login node as a job can take days to finish. A simple tcsh script (.csh) to submit job is attached in '/example/'. There are three steps to run a new simulations:

4.1 Generate initial configuration. This step is done in serial. Replace the bold line by the real path to *initconfig* on your system.
> /**path_to_executive_file_initconfig**/initconfig

For example: If you save the package to '/home/user/Parallel-DMD-PRIME20' then the path to executable file will be '/home/user/Parallel-DMD-PRIME20/src/'. Your submission script will look like:

> **/home/user/Parallel-DMD-PRIME20/src**/initconfig

Before a DMD simulation, the system will be heated to a high temperature and then be slowly annealed to the desired temperature. This step is to make sure that all peptide chains are denatured and that the DMD simulation starts with all random coils. The submission command for annealing process is as follow.
> foreach i (`seq 1 number_of_temperatures_use_for_annealing`)
> 
> mpirun /**path_to_executive_file_DMDPRIME20**/DMDPRIME20 < inputs/annealtemp_$i > outputs/out_annealtemp_$i
> 
> end

If using default temperature, *number_of_temperatures_use_for_annealing = 9*

If using user-defined temperature:



The simulation temperature and numbers of collisions are defined by users. Larger system will need longer simulation times. The frequency to write ouput is defined in *input.txt*. This value is called a round of simulation. The total number of collisions will depend on number of simulation rounds which is defined in submission script.
> foreach i (`seq 1 number_of_simulation_round`)
> 
> mpirun /**path_to_executive_file_DMDPRIME20**/DMDPRIME20 < inputs/simtemp_$i > outputs/out_simtemp_$i
> 
> end

## Developing Status
The software is being developed and updated. An result analysis package is being developed.

## References:
[1] Wang, Y., Shao, Q., and Hall, C. K. N-terminal Prion Protein Peptides (PrP(120–144)) Form Parallel In-register β-Sheets via Multiple Nucleation-dependent Pathways. Journal of Biological Chemistry. Vol. 292, Issue 50. (2016). https://doi.org/10.1074/jbc.M116.744573

[2] Wang, Y., Shao, Q., and Hall, C. K. Seeding and cross-seeding fibrillation of N-terminal prion protein peptides PrP(120–144). Protein Science (2018). https://doi.org/10.1002/pro.3421

[3] Cheon, M., Chang, I., and Hall, C. K. Influence of temperature on formation of perfect tau fragment fibrils using PRIME20/DMD simulations. Protein Science (2012). https://doi.org/10.1002/pro.2141

[4] Cheon, M., Chang, I., and Hall, C. K. Spontaneous Formation of Twisted Ab16-22 Fibrils in Large-Scale Molecular-Dynamics Simulations. Biophysical Journal. Vol. 101, 2493-2501 (2011).  https://doi.org/10.1016%2Fj.bpj.2011.08.042

[5] Samuel J. Bunce et al., Molecular insights into the surface-catalyzed secondary nucleation of amyloid-β40 (Aβ40) by the peptide fragment Aβ16–22.Sci. Adv.5,eaav8216(2019).DOI:10.1126/sciadv.aav8216

[6] Yiming Wang et al., Thermodynamic phase diagram of amyloid-β (16–22) peptide. Proceedings of the National Academy of Sciences. Vol. 116, Issue 6, 2091-2096 (2019). doi:10.1073/pnas.1819592116

[7] Wang, Y., Latshaw, D. C., and Hall, C. K. Aggregation of Aβ(17–36) in the Presence of Naturally Occurring Phenolic Inhibitors Using Coarse-Grained Simulations. Journal of Molecular Biology. Volume 429, Issue 24, 3893-3908 (2017). https://doi.org/10.1016/j.jmb.2017.10.006.

[8] Xingqing Xiao and others, Sequence patterns and signatures: Computational and experimental discovery of amyloid-forming peptides, PNAS Nexus, Volume 1, Issue 5, November 2022, pgac263, https://doi.org/10.1093/pnasnexus/pgac263
