function [Ind]=DrawCircle(Node,Element)

%ChooseCircle Lets you to define a circular inhomogeneity for 2D EIT simulations
% Function [Ind]=ChooseCircle(Node,Element);
% lets you to define a circular inhomogeneity for 2D EIT
% simulations. Choose the center with the left mouse button and define the radius with
% the right button.
%
% INPUT
%
% Node = nodal data structure
% Element = element data structure
%
% OUTPUT
%
% Ind = Indices of the chosen elements


% M. Vauhkonen 17.4.2000
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


Nod=reshape([Node.Coordinate]',2,max(size(Node)))';
Topol=reshape([Element.Topology]',size(Element(1).Topology,2),max(size(Element)))';
clf,plcigrid(Nod,Topol,[],'r'),axis equal;
NTopol=max(size(Topol));

% Choose the center and the radius
[x,y]=getline(gcf);
r=norm([x(1),y(1)]-[x(2),y(2)]);
Ind = sparse(zeros(NTopol,1));

for ij=1:NTopol
 nodes=Nod(Topol(ij,:),:);
 center = mean(nodes);
    if norm(center(:)-[x(1);y(1)])<r
       Ind(ij) = 1;
    end
end
Ind=find(Ind);
 for ii=1:max(size(Ind))
  Hii=Nod(Topol(Ind(ii),:),:);
  patch(Hii(:,1),Hii(:,2),1);
 end




