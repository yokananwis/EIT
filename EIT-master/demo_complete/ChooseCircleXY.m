function [Ind] = ChooseCircleXY(Node,Element,point1,point2)
%CHOOSECIRCLEXY is the same as ChooseCircle, but the user provides the
% endpoints for the lines defining the circle radius, in the form point=[x,y]. 

Nod=reshape([Node.Coordinate]',2,max(size(Node)))';
Topol=reshape([Element.Topology]',size(Element(1).Topology,2),max(size(Element)))';
clf,plcigrid(Nod,Topol,[],'r'),axis equal;
NTopol=max(size(Topol));


r=norm(point1-point2);
Ind = sparse(zeros(NTopol,1));

for ij=1:NTopol
 nodes=Nod(Topol(ij,:),:);
 center = mean(nodes);
    if norm(center(:)-[point1(1);point1(2)])<r
       Ind(ij) = 1;
    end
end
Ind=find(Ind);
 for ii=1:max(size(Ind))
  Hii=Nod(Topol(Ind(ii),:),:);
  patch(Hii(:,1),Hii(:,2),1);
 end

end

