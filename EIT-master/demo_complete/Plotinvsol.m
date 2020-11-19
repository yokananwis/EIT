function []=Plotinvsol(rho,g,H)

%Plotinvsol Plots the solution of 2D EIT problem with piecewise constant basis
% Function []=Plotinvsol(rho,g,H);
% plots the solution of 2D EIT problem.
%
% INPUT
%
% rho = the solution, constant on each element
% g = node coordinates (x,y) (N x 2 matrix)
% H = connectivity matrix 

% J. Kaipio, 11.4.1994. Modified by J. Virkkala 19.9.1994
% Modified for EIDORS by M. Vauhkonen 13.8.1999
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


[nH,~]=size(H);

axis('off'),axis('square')
set(gcf,'defaultpatchedgecolor','none');

for ii=1:nH
  Hii=g(H(ii,:),:);
   patch(Hii(:,1),Hii(:,2),rho(ii));		
end









