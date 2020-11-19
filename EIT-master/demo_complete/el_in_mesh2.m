function Ind = el_in_mesh2(H1,g1,H2,g2)

%el_in_mesh2 Makes the element connectivity data matrix between the two meshes 
% Function Ind = el_in_mesh2(H1,g1,H2,g2) 
% makes the element connectivity data matrix between the two meshes  defined by 
% H1, g1 and H2, g2. 
%
% INPUT
%
% H1 = connectivity of the coarse mesh
% g1 = coordinates of the nodes of the coarse mesh
% H2 = connectivity of the dense mesh
% g2 = coordinates of the nodes of the dense mesh
%
% OUTPUT
%
% Ind = element connectivity data matrix
 
% L.M. Heikkinen 5.5.1999
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: laheikki@venda.uku.fi

S1 = size(H2,1);
S2 = size(H1,1);
Ind = zeros(S1,S2);
for i = 1:S1
  g = g2(H2(i,:),:);
  g11 = (1/3)*sum(g,1);
  T = tsearch(g1(:,1),g1(:,2),H1,g11(1,1),g11(1,2));
  Ind(i,T) = 1;
end
Ind = sparse(Ind);


