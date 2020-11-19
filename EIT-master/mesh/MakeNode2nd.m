function [Node]=MakeNode2nd(Element,Nodelist,g);

%MakeNode2nd Computes the Node data structure for a 2D mesh having quadratic triangles
% Function [Node]=MakeNode2nd(Element,g);
% computes the Node data for MeshData.
% Node is a structure including all the nodal coordinates and
% for each node there is information to which nodes (NodeConnection) 
% and elements (ElementConnection) the node is
% connected.  
%
% INPUT
%
% Element = Element structure, see MakeElement.m
% Nodelist = list of all the faces, see MakeElement.m
% g = coordinates of the nodes
%
% OUTPUT
%
% Node = node structure

% M. Vauhkonen, University of Kuopio, Finland, 11.8.1999
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

[rg,cg]=size(g);
msE=max(size(Element));
for ii=1:rg
 ElementConnection=[];
 Node(ii).Coordinate=[g(ii,:)];
  for jj=1:msE
   if find(Element(jj).Topology==ii)
    ElementConnection=[ElementConnection,jj];
   end
  end
  Node(ii).ElementConnection=ElementConnection;
 Nc=[];
 [I,J]=find(Nodelist==ii);
 apu1=Nodelist(I(find(J==1)),:); %Corner node
 apu2=Nodelist(I(find(J==2)),:); %Middle node
 apu3=Nodelist(I(find(J==3)),:); %Corner Node
 apu=[apu1(:,2);apu2(:,1);apu2(:,3);apu3(:,2)];
 Nc=unique(apu);
  Node(ii).NodeConnection=Nc;
end
  



