%MakeMeshData2nd A script that builds two circular meshes, quadratic triangles,  for 2D EIT Package demos 
% MakeMeshData2nd is a script that builds two circular  
% meshes for 2D EIT Package. Denser mesh (mesh 2) is for forward computations and
% a coarse mesh (mesh 1) for inverse computations. See also
% meshgen_eit2d.m, MakeElement2nd and MakeNode. 
% QMG is used for mesh generation, see http://www.cs.cornell.edu/Info/People/vavasis/qmg-home.html.


% M. Vauhkonen, 11.8.1999
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


S=2.5; % Length of the electrode.
N=16;  % Number of the electrodes.
r=14;  % Radius of the circle.
style='s'; % 's' for the structured mesh and 'u' for the unstructured mesh

[H1,g1,H2,g2,E1,E2,Ind2,Indb1,Indb2] = meshgen_eit2d(S,N,r,style);
g=g1;H=H1;                                  %Old data for plotting etc.
[Element,Nodelist]=MakeElement(H1,Indb1,E1); %-----------"------------
[Node]=MakeNode(Element,Nodelist,g);        %-----------"----------  

[g1,H1,Indb1]=addnodes(g1,H1,Indb1);
[Element1,Nodelist1]=MakeElement2nd(H1,Indb1,E1);
[Node1]=MakeNode2nd(Element1,Nodelist1,g1);




