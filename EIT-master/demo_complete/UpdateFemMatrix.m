function [A]=UpdateFemMatrix(Agrad,Kb,M,S,sigma);

%UpdateFemMatrix Assembles the system matrix for EIT
% Function [A]=UpdateFemMatrix(Agrad,Kb,M,S,sigma);
% updates the system matrix A given a new conductivity
% vector sigma.
%
% INPUT
%
% Agrad = the gradient part of the system matrix (see FemMatrix.m)
% Kb,M and S = other blocks of the system matrix (see FemMatrix.m)
% sigma = conductivity (or admittivity) vector
%
% OUTPUT
%
% A = the updated system matrix 

% M. Vauhkonen 23.3.2000, University of Kuopio, Department of Applied Physics,
% PO Box 1627, FIN-70211, Kuopio, Finland, email:Marko.Vauhkonen@uku.fi

rAgrad=sqrt(size(Agrad,1));
A=[reshape(Agrad*sigma,rAgrad,rAgrad)+Kb,M;M.',S];





