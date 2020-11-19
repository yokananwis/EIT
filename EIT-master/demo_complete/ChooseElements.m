function Ind=ChooseElements(Node,Element,E,N)

%ChooseElements Lets you to choose an inhomogeneity (a set of elements) for 2D EIT simulations
% Function Ind=ChooseElemenst(Node,Element,E,N);
% lets you choose N elements from a mesh
% defined by Element. If N is not given, as many
% elements as you want are chosen. Ind are the indices
% to the Element topolgy.
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
  % I don't think this actually works because it's ignoring the
  % pre-determined topology:
  Ind(ii) = pointLocation(d,q);
  %Ind(ii)=tsearch(Nod(:,1),Nod(:,2),Topol,x,y);
  Hii=Nod(Topol(Ind(ii),:),:);
  patch(Hii(:,1),Hii(:,2),1);            
 end
else
  Ind=[];
  [x,y]=ginput;
  d = delaunayTriangulation(Nod(:,1),Nod(:,2));
  q = [x,y];
  Ind = [Ind;pointLocation(d,q)];
  %Ind=[Ind;tsearch(Nod(:,1),Nod(:,2),Topol,x,y)];
 for ii=1:max(size(Ind))
  Hii=Nod(Topol(Ind(ii),:),:);
  patch(Hii(:,1),Hii(:,2),1);           
 end
end



