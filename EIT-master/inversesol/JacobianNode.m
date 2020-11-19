function [J]=JacobianNode(Node,Element,Agrad,U,U0,rho,style);

%JacobianNode Computes the Jacobian for 2D EIT when node basis is used
% Function [J]=JacobianNode(g,H,Agrad,U,U0,rho,style);
% computes the Jacobian for 2D EIT when node basis
% is used.
%
% INPUT
%
% Node = nodal data structure
% Element = element data structure
% Agrad = \int_{Element(ii) \nabla\phi_i\cdot\nabla\phi_j
% U = voltages corresponding to the injected currents
% U0 = voltages corresponding to the measurement field
% rho = resistivity or admittivity vector
% style = either 'real' for reconstructing resistivity or 'comp'
% for reconstructin admittivity
%
% OUTPUT
%
% J = Jacobian

% M. Vauhkonen 13.8.1999,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


NNode=max(size(Node));
NElement=max(size(Element));



if nargin<3

mE=max(size(Element(1).Topology));
if mE==3
 Agrad=sparse(NNode^2,NNode);
 H=reshape([Element.Topology],3,NElement)';
 mH=max(max(H)); 
else
 H=reshape([Element.Topology],6,NElement)';
 mH=max(max(H(:,1:2:6))); 
 Agrad=sparse(NNode^2,mH);
 clear H
end
g=reshape([Node.Coordinate],2,NNode)';      %Nodes

for jj=1:mH
Aa=sparse(NNode,NNode);
El=Node(jj).ElementConnection;
 for ii=1:max(size(El))
   ind=Element(El(ii)).Topology; % Indices of the element
   gg=g(ind,:);
   if max(size(gg))==3
    indsig=ind;
    I=find(jj==indsig);
    anis=grinprodgausnode(gg,I);
    Aa(ind,ind)=Aa(ind,ind)+anis;
   else
     indsig=ind(1:2:6);
     I=find(jj==indsig);
     anis=grinprodgaus2ndnode(gg,I);
     Aa(ind,ind)=Aa(ind,ind)+anis;         
   end  
 end  
 Agrad(:,jj)=Aa(:);
end
J=Agrad;
else
if style=='comp'
  J=zeros(size(U,2)*size(U0,2),2*size(Agrad,2));
 for ii=1:size(Agrad,2);
    Agrad1=reshape(Agrad(:,ii),NNode,NNode);
    AgU1=Agrad1*U(1:NNode,:);
    AgU2=Agrad1*U(NNode+1:2*NNode,:);
    JJ=-U0.'*[AgU1;AgU2];
    JJ=JJ(:);
    J(:,ii)=JJ;
    JJ=-U0.'*[-AgU2;AgU1];
    J(:,ii+size(Agrad,2))=JJ(:);
 end
elseif style=='real'
  J=zeros(size(U,2)*size(U0,2),size(Agrad,2));
   for ii=1:size(Agrad,2);
    JJ=-U0.'*reshape(Agrad(:,ii),NNode,NNode)*U;
    JJ=JJ(:);
    J(:,ii)=JJ;
   end
 end
end








