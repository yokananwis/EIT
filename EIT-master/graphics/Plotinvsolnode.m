function [fh,fp]=Plotinvsolnode(sol,g,H);

%Plotinvsolnode Plots the solution of 2D EIT problem in linear basis 
% Function []=Plotinvsolnode(sol,g,H);
% plots the solution of 2D EIT problem. Conductivity in linear basis.
%
% INPUT
%
% sol = the solution, values in the nodes
% g = node coordinates (x,y) (N x 2 matrix)
% H = connectivity matrix

% J. Kaipio, 11.4.1994. Modified by J. Virkkala 19.9.1994
% Modified for EIDORS by M. Vauhkonen 13.8.1999
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi



[solN,solM]=size(sol);
[gN,gM]=size(g);
[HN,HM]=size(H);

for ii=1:solM
  view(2);
  hold on
  for ij=1:HN
    X=g(H(ij,:),1);
    Y=g(H(ij,:),2);
    Z=sol(H(ij,:),ii);
    fp=patch(X,Y,Z,Z);
    set(fp,'edgecolor','none');
    hold on
  end    
end
axis('square')









