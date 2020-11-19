function indbound = findboundary(g,H)

%findboundary Indices of the boundary nodes
% Function indbound = findboundary(g,H);
% finds the indices of the boundary nodes.
%
% INPUT
% 
% g = node coordinates
% H = connectivity
%
% OUTPUT
% 
% indbound = indices of the boundary nodes

% L.M. Heikkinen 2.6.1999 
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: laheikki@venda.uku.fi


indbound = [];
for i = 1:size(g,1);
 I11 = find(H == i);
  if size(I11,1) < 4
     indbound = [indbound;i];
  end
end                                      
