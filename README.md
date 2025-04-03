Visualizing Influenza Virus Evolutionary Trajectories in Multidimensional Landscapes

## Problem Statement
The primary challenge is to simplify and visualize complex sequence surveillance data to assist in leadership decision-making and inform public health policies effectively.

Solution Approach
The project employs fitness landscape modeling to depict complex virological data relationships. A 3D modeled surface, created by overlaying a two-dimensional plane of continuous data distributions with an independent "fitness" variable, is sliced into 2D plots along a summary path using locally-weighted regression. This method produces visualizations that map the fitness peaks and troughs of the data, thereby illustrating evolutionary "fitness" of various viral strains.

## Project Structure


bin/: Contains executable scripts and binaries.

sql/: SQL scripts for database management and queries.

FitnessLandscape: Bash wrapper
