function OutSpMat = SparseCrush(SpMat)

%SparseCrush  Minimise memory requirement of sparse matrix
% Function OutSpMat = SparseCrush(SpMat);
% Minimises memory requirement of a sparse matrix.
%
% INPUT 
% SpMat = sparse matrix with unefficient use of memory
%
% OUTPUT
%
% OutSpMat = same as SpMat but uses less memory (in certain cases)    

% Dr Richard F Booth 2000,
% Dept. of Maths,  UMIST, 
% PO Box 88, Manchester M60 1QD, email: richard.booth@umist.ac.uk

[ i, j, s ] = find( SpMat ) ;
[ m, n ] = size( SpMat ) ;
OutSpMat = sparse( i , j , s , m , n ) ;
