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
