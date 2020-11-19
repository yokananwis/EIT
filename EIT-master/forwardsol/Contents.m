% Functions for solving the EIT forward problem in 2D
%
% FEM integrals for the linear triangular elements
%
%  grinprodgaus     - Computes the integral of the product of the gradients in 2D FEM
%  bound1           - Computes the boundary integral of one linear basis function in 2D FEM
%  bound2           - Computes the boundary integral of the product of two linear basis functions in 2D FEM
%
% FEM integrals for the isoparametric quadratic triangular elements
%
%  grinprodgausquad - Computes the integral of the product of the gradients in 2D FEM for quadratic
%                     isoparametric triangular elements
%  grinprodgausnode - Computes the gradient part of the system matrix
%                     in the linear basis case. The conductivity  is also in the linear basis.
%  grinprodgausnodequad  - Computes the gradient part of the system matrix
%                         in the quadratic basis case. The conductivity  is in the linear basis.
%  boundquad1       - Computes the boundary integral of one quadratic basis function 
%  boundquad2       - Computes the boundary integral of the product of two quadratic basis functions 
%  electrlen        - Computes the length of the electrode in the case of isoparametric (quadratic) triangular element
%
% Functions for forward computations
%
%  Current          - Supplies some current patterns for 2D EIT
%  FemMatrix        - Computes the blocks of the system matrix for 2D EIT with linear and quadratic basis 
%  FemMatrixNode    - Computes the blocks of the system matrix for 2D EIT with linear and quadratic basis. 
%                     Conductivity in linear basis.
%  ForwardSolution  - Solves the potential distribution in 2D EIT
%  RemoveVolt       - Removes all the voltages measured on the current carrying electrodes
%  UpdateFemMatrix  - Assembles the system matrix for EIT
%
% Sparse trick
%
%  SparseCrush      -  Minimise memory requirement of sparse matrix







