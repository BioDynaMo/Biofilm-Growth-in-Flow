# Biofilm Growth in Flow

This project provides a three-dimensional ABM of a biofilm growing under constant fluid flow that captures the phenomena of erosion, sloughing and dispersion throughout the biofilm life cycle. A multi-scale approach was implemented to include the long timescale of growth and EPS production and the short timescale of bacterial collisions and detachment. These simulations have been validated both qualitatively and quantitatively with experimental data.

This project is currently being updated to reduce bloat in the code. Currently this version of the project only runs on linux due to the implementation of OpenFoam coupling with BioDynaMo. We recommend running this project on a high-performance computer as the particle-laden flow portion of the simulation is computationally demanding.

## Prerequisites

Install BioDynaMo and OpenFoam10. It is recommended to install BioDyanMo first as installing openfoam10 automatically installs ParaviewOpenFoam10 which causes a BioDynaMo to fail installation. Once these two softwares have been installed, you can run the project.

## Commands to run the project
Once BioDynaMo and OpenFoam are installed, open a separate terminal and download the project in a place of your choosing. Enter the directory and run the following commands.

Run ```source BioDynaMo``` to source the biodynamo code and then run ```bdm build```.

Run ```source bashrc``` to source OpenFoam without sourcing paraFoam. Then run ```wmake updateGrid``` to build the updateGrid function used to create the fluid flow field.

Then run ```bdm run``` to run the simulation. Currently a simulation will take 2 - 4 days to run so it is recommended to either run the simulation on a separate computer or partition BioDynaMo to use only a portion of the computer.


## File breakdown

OpenFoam Files:  
```0``` - stores the OpenFoam field values for fluid velocity and pressure  
```constant``` - stores the OpenFoam FVM mesh where the stokes equation is calculated  
```system``` - Extra files required by OpenFoam to solve the stokes equation  

Coupling files:  
```cfd_grid.csv``` - Contains (x, y, z) coordinates for each voxel center in the FVM mesh in BioDynaMo  
```cfdOutput.csv``` - Contains fluid velocity values corresponding to each voxel in the FVM mesh in BioDynaMo  
```updateGrid``` - Custom OpenFoam function which reads in cfd_grid.csv and outputs cfdOutput.csv  

Files needed to run the project:  
```CMakeList.txt``` - Needed to build the project source code  
```test``` - Needed to build the project source code  
```bashrc``` - Used to source OpenFoam without paraFoam  
```allClean``` - Deletes all output files created during the simulation  
```bdm.json``` - Used to generate visualisation files for paraview  

Project source code:  
```src``` - Contains the project source code  
```initial_chunk.csv``` - Information for initial agents to begin the simulation. Each row is one agent and the columns are labelled as positions (x, y, z), Volume of cell, Volume of EPS shell, cell density.  

## Source Code breakdown

All source code files are found within ```src```. All agents are defined with the class ```MyCell```, cells are defined with the class ```HET``` and EPS particles are defined with the class ```EPS```.  

```biofilm_growth.h``` and ```biofilm_growth.cc``` - Main project file where the agents are initialised, the scheduler is set up and the simulation runs.  
```sim_param.h``` - Contains custom parameters for simulation.   
```cell_class/mycell.h``` - Base class for all agents, giving all agents a velocity, net force and other variables required for calculating forces. Contains functions for calculating agent-agent force, agent-surface force and drag force, and updating velocity.  
```cell_class/HET.h``` - Extension of ```MyCell``` class to define cells. Includes definitions and functions which separate volume of cell and volume of EPS shell.  
```cell_class/EPS.h``` - Extension of ```MyCell``` class to define EPS particles. Only used for visualisation purposes.  
```biology_functions/Growth.h``` - Growth and division behaviour given to cells.  
```biology_functions/MakeEPS.h``` - EPS production behaviour given to cells.  
```biology_functions/sensing.h``` - Quorum sensing behaviour given to cells if switched on.  
```RelaxAgents.h``` - Operation to relax the biofilm before initiating the particle-laden flow. Loop which calculates net force, applies net force and checks if net force on all agents is below 0.5 nN.  
```PhysicalSteadyState.h``` - Operation to run particle-laden flow until steady-state after biofilm has relaxed. Loop which calculates net force, applies net force and checks if all agents are either static or outside the simulation.  
```Boundary.h``` - Contains function which calculates force from boundaries (this does not include maximum x boundary as it is considered open, and it does not include minimum z boundary as this is the surface, where the agent-surface forces are calculated in ```cell_class/mycell.h```).  
```MyMath.h``` - Contains functions to calculate mod and dot product of two vectors.
```CFDgrid.h``` and ```CFDgrid.cc``` - Class defining the FVM mesh which contains the fluid flow field in BioDynaMo. Initialising the grid produces the ```cfd_grid.csv``` file and reads in the ```cfdOutput.csv``` file. Class includes retrieving the fluid flow velocity at given position and and calculating the solid volume fraction field.  

## Units

In order to appropriately generate the paraview files, all parameters are converted into custom units which are related to SI units using the following conversion:  
Length: 1 unit = 1 um = 1e-6 meters  
Mass: 1 unit = 1e-18 kilograms  
Time: 1 unit = 550 seconds for biological timestep, 1 unit = 1e-8 seconds for physical timestep.
