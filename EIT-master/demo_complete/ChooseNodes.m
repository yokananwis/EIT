function Ind=ChooseNodes(Node,Element,E,N)

%ChooseNodes Lets you to choose an inhomogeneity (a set of nodes) for 2D EIT simulations 
% Function Ind=ChooseNodes(Node,Element,E,N);
% lets you choose N elements from a mesh
% defined by Node and Element. If N is not given, as many
% elements as you want are chosen. 
%
% INPUT
%
% Node = node information
% Element = element information
% E = electrode information (elements under the electrodes, optional, use [] if not defined)
% N = the number of elements to be chosen (optional)
%
% OUTPUT
%
% Ind = Indices of the chosen elements

% M. Vauhkonen 17.8.1999,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


Nod=reshape([Node.Coordinate]',2,max(size(Node)))';
Topol=reshape([Element.Topology]',size(Element(1).Topology,2),max(size(Element)))';
clf,plcigrid(Nod,Topol,E,'r'),axis equal;
hold on
if nargin==4
 Ind=zeros(N,1);
 for ii=1:N
  [x,y]=ginput(1);
  d = delaunayTriangulation(Nod(:,1),Nod(:,2));
  q = [x,y];
  Ind(ii) = nearestNeighbor(d,q);
  %Ind(ii)=dsearch(Nod(:,1),Nod(:,2),Topol,x,y);
  plot(Nod(Ind(ii),1),Nod(Ind(ii),2),'*')
 end
else
 Ind=[]; 
 [x,y]=ginput;
 d = delaunayTriangulation(Nod(:,1),Nod(:,2));
 q = [x,y];
 Ind = nearestNeighbor(d,q);
 %Ind=[Ind;dsearch(Nod(:,1),Nod(:,2),Topol,x,y)];
 for ii=1:max(size(Ind))
   plot(Nod(Ind(ii),1),Nod(Ind(ii),2),'*')
 end
end

%Ind=dsearch(Nod(:,1),Nod(:,2),Topol,x,y);

