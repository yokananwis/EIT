% Functions for mesh handling
%
% Plotting
%
%  plcigrid        - Plots a 2D FEM mesh
%
% Making inhomogeneities
%
%  ChooseCircle    - Lets you to define a "circular" inhomogeneity for 2D EIT simulations
%  ChooseElements  - Lets you to choose an inhomogeneity (a set of elements) for 2D EIT simulations
%  ChooseNodes     - Lets you to choose an inhomogeneity (a set of nodes) for 2D EIT simulations 
%
% Making meshes, general
%
%  meshgen_eit2d   - Makes two circular meshes. You need the QMG mesh generator
%  findboundary    - Indices of the boundary nodes
%  cirgrid_eit     - Makes circular meshes for 2D EIT
%  el_in_mesh2     - Makes the element connectivity data matrix between the two meshes 
%  RefineMesh      - Refines a coarse mesh by adding new nodes in the triangle edges 
%
% Making meshes, linear triangles
%
%  MakeMeshData    - A script that builds two circular meshes for 2D EIT Package demos
%  MakeElement     - Computes the Element data structure for a 2D mesh
%  MakeNode        - Computes the Node data structure for a 2D mesh
%
% Making meshes, quadratic triangles
%
%  MakeMeshData2nd - A script that builds two circular meshes, quadratic triangles,  for 2D EIT Package demos 
%  MakeElement2nd  - Computes the Element data structure for a 2D quadratic triangular mesh
%  MakeNode2nd     - Computes the Node data structure for a 2D mesh having quadratic triangles
%  addnodes        - Adds nodes in the middle of each face of the triangular element






