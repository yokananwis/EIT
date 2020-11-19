%MakeMeshData A script that builds two circular meshes for 2D EIT Package demos 
% MakeMeshData is a script that builds two circular  
% meshes for 2D EIT Package demos. Denser mesh (mesh 2) is for forward computations and
% a coarse mesh (mesh 1) for inverse computations. See also 
% meshgen_eit2d.m , MakeElement and  MakeNode. 
% QMG is used for mesh generation, see http://www.cs.cornell.edu/Info/People/vavasis/qmg-home.html. 

% L.M. Heikkinen and M. Vauhkonen, 11.8.1999
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

S=2.5; % Length of the electrode.
N=16;  % Number of the electrodes.
r=14;  % Radius of the circle.
style='s'; % 's' for the structured mesh and 'u' for the unstructured mesh

[H1,g1,H2,g2,E1,E2,Ind2,Indb1,Indb2] = meshgen_eit2d(S,N,r,style);
[Element1,Nodelist1]=MakeElement(H1,Indb1,E1);
[Node1]=MakeNode(Element1,Nodelist1,g1);
[Element2,Nodelist2]=MakeElement(H2,Indb2,E2);
[Node2]=MakeNode(Element2,Nodelist2,g2);


%% If you do not have access to QMG, try next one

r=[14,12,10,8,5,3,0];
N=[32,26,20,16,12,8,1];
%N=[128,100,70,45,20,10,1];
eI=[1,1]; 
[g1,gp,H1,E1]=cirgrid_eit(r,N,eI);
[g2,H2,Ind2] = RefineMesh(g1,H1);

E2=[];
for ii=1:size(E1,1)
 Eii=E1(ii,:);
 E2ii=[];
  for jj=1:size(Eii,2)
   Ind=find(Ind2(:,Eii(jj)));
   E2ii=[E2ii,Ind(1:2)']; 
  end
 E2=[E2;E2ii];
end 

Indb1=findboundary(g1,H1);
Indb2=findboundary(g2,H2);

[Element1,Nodelist1]=MakeElement(H1,Indb1,E1);
[Node1]=MakeNode(Element1,Nodelist1,g1);
[Element2,Nodelist2]=MakeElement(H2,Indb2,E2);
[Node2]=MakeNode(Element2,Nodelist2,g2);






