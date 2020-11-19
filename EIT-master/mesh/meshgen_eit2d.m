function [H1,g1,H2,g2,E1,E2,Ind2,Indb1,Indb2] = meshgen_eit2d(S,N,r,style)

%meshgen_eit2d Makes two circular meshes. You need the QMG mesh generator
% Function [H1,g1,H2,g2,E1,E2,Ind2,Indb1,Indb2] = meshgen_eit2d(S,N,r,style)
% makes two  circular meshes, a coarse (submesh of the dense mesh)
% and a dense meshes, for 2D EIT given the size and the number of the electrodes,
% radius of the tank and the style of the mesh.
% See also  coord_for_elec.m, findboundary.m and el_in_mesh2.m
%
% INPUT
%
% S = size of the electrode
% N = number of the eletrodes
% r = radius of the tank
% style = 's' for the structured mesh and 'u' for the unstructured mesh
%
% OUTPUT
%
% H1 = connectivity of the coarse mesh
% g1 = coordinates of the nodes of the coarse mesh
% H2 = connectivity of the dense mesh 
% g2 = coordinates of the nodes of the dense mesh
% E1 = electrode information for the coarse mesh
% E2 = electrode information for the dense mesh
% Ind2 = element connectivity data matrix between the two meshes
% Indb1 = boundary node indices for the coarse mesh
% Indb2 = boundary node indices for the dense mesh

% L.M. Heikkinen 10.5.1999 & 1.6.1999
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: laheikki@venda.uku.fi

[x,y] = coord_for_elec(S,N,r);
vtx = [x y];

[p,l,s,c] = gmbasics;
point = gmembed(gmembed(p));  

brep = gm_poly(vtx,0);
br = gmpolygon(32);
br = gmapply(gmdilate(r-1,r-1),br);
brep1 = ed_addregion(brep,br);

if style == 'u'
 mesh = ed_meshgen(brep1,3);
end

if style == 's'
 i = 10;
 j = 2;
 while j <= (r - 2),
   br = gmpolygon(i);
   br1 = gmapply(gmdilate(j,j),br);
   brep = ed_addregion(brep,br1);
   i = i + 6;
   j = j + 2; 
 end
 br = gmpolygon(size(vtx,1));
 br1 = gmapply(gmdilate(r-1,r-1),br);
 brep = ed_addregion(brep,br1);
 brep = ed_addregion(brep,point);
 mesh = ed_meshgen(brep,3);
end

mesh2 = gmrefine(mesh);

g1 = gmget(mesh,'vertices');
H1 = gmget(mesh,'simplices')+1;

g2 = gmget(mesh2,'vertices');
H2 = gmget(mesh2,'simplices')+1;

%% Elements under the electrodes 

Indh = 1:size(x,1);
Indh = reshape(Indh,4,N)';
Ind2 = feval('el_in_mesh2',H1,g1,H2,g2);
E1 = [];
E2 = [];
E = 0;
for i = 1:N  
 for j = 1:3
    E = find(H1(:,1) == Indh(i,j) & H1(:,2) == Indh(i,j+1));
    E1 = [E1;E];
    Eh = find(Ind2(:,E)); 
    E2 = [E2;Eh(1:2)]; 
 end
end

E1 = reshape(E1,3,N);
E2 = reshape(E2,6,N);
E1 = E1';
E2 = E2';

Indb1 = feval('findboundary',g1,H1);
Indb2 = feval('findboundary',g2,H2);



