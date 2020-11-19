function [J]=Jacobian(Node,Element,Agrad,U,U0,rho,style);

%Jacobian Computes the Jacobian for 2D EIT when elementwise basis is used
% Function [J]=Jacobian(g,H,Agrad,U,U0,rho,style);
% computes the Jacobian for 2D EIT when elementwise basis
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
Agrad=sparse(NNode^2,NElement); % Gradients of the basis functions integrated over 
                       % each element. 

 for ii=1:NElement
  ind=Element(ii).Topology;
  gg=reshape([Node(ind).Coordinate],2,max(size(ind)))'; % A 3x2 or 6x2 matrix of triangle vertices in (x,y) coord.
  if max(size(ind))==3
   anis=grinprodgaus(gg,1);
  else
   anis=grinprodgausquad(gg,1);
  end
  Aa=sparse(NNode,NNode);
  Aa(ind,ind)=anis;
  Agrad(:,ii)=Aa(:);    
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
  JJ=U0.'*1/rho(ii)^2*reshape(Agrad(:,ii),NNode,NNode)*U;
  JJ=JJ(:);
  J(:,ii)=JJ;
 end
end
end


