# EIT
Electrical Impedance Tomography Reconstruction Algorithm Based on Finite Element Model Using the Eidors Library in Octave/Matlab
Eidors website: http://eidors3d.sourceforge.net/

How to Use EIDORS:
1. Unzip the Eidors software into the designated directory.
2. Start Matlab and run the command run /path/to/eidors/startup.m.
3. Try running the programs eidors_demo or eidors_demo2 in the EIT-master folder.

# Simulation
This folder contains simulation programs for reconstruction using the Newton-Raphson method (Regularization, Pseudoinverse, Singular Value Decomposition) and the Landweber Iteration method.
To try the simulations, copy the simulation files (except ChooseCircle.m) to the EIT-master\demo_complete directory.
The ChooseCircle.m file has been updated to be compatible with Octave, which does not support defining the position and size of circles.

# Reconstruction Program
This folder contains the EIT reconstruction program using real boundary voltage data from an Excel file.
To try this EIT reconstruction program, copy the program files and the Excel data files to the EIT-master\demo_complete directory.

# Mesh Generator
The mesh generator uses Netgen 6.1, which requires the following installations:
1. Python 3.7
2. Microsoft visual studio 2017 redistributable
3. EIDORS
It is recommended to install Netgen in the EIDORS directory.
Documentation for the functions in Netgen can be found at: http://eidors3d.sourceforge.net/doc/index.html?eidors/meshing/netgen/call_netgen.html.
Examples of mesh creation are available in the Netgen folder.
Netgen website: https://ngsolve.org/
