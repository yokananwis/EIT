function []=plcigrid(g,H,E,co)

%plcigrid Plots a 2D FEM mesh  
% Function []=plcigrid(g,H,E,co);
% plots a given 2D FEM mesh.
%
% INPUT
%
% g = node coordinate matrix
% H = topology 
% E = elements under the electrodes (can also be empty)
% co = colour of the mesh

% J. Kaipio, 11.4.1994. 
% Modified by M. Vauhkonen 1999,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


nH=max(size(H));

axis('xy'),axis('square')
hold on


for ii=1:nH
  Hii=g(H(ii,:),:);
  Hii=[Hii;Hii(1,:)];
  if ~isempty(E)
    hHii=plot(Hii(:,1),Hii(:,2),co);
    [pp,cc]=find(ii==E);
    if ~isempty(find(ii==E,1)),
     hHii=fill(Hii(:,1),Hii(:,2),pp);
    end
  else
    hHii=plot(Hii(:,1),Hii(:,2),co);
  end
end












