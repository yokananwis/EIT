function [Agrad,Kb,M,S,C]=FemMatrixNode(Node,Element,z);

%FemMatrixNode Computes the blocks of the system matrix for 2D EIT with linear and quadratic basis. Conductivity in linear basis.
% Function [Agrad,Kb,M,S,C]=FemMatrixNode(Node,Element,z);
% computes the matrices needed in the finite element
% approximation of the 2D EIT forward problem. Conductivity in linear basis.
%
% INPUT
%
% Node = nodal data structure
% Element = element data structure
% z = a vector of (complex) contact impedances
%
% OUTPUT
%
% Agrad = the gradient part of the system matrix
% Kb,M and S = other blocks of the system matrix
% C = voltage reference matrix

% M. Vauhkonen 11.5.1994, modified from the version of J. Kaipio
% 25.4.1994. Modified 5.9.1994 by M. Vauhkonen for EIT.
% Modified for EIDORS by M. Vauhkonen 11.5.2000
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


Nel=max(size(z));                           %The number of electrodes.
NNode=max(size(Node));                      %The number of nodes
NElement=max(size(Element));                %The number of elements
M=sparse(NNode,Nel);
Kb=sparse(NNode,NNode);
s=zeros(Nel,1);
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
    anis=grinprodgausnodequad(gg,I);
    Aa(ind,ind)=Aa(ind,ind)+anis;
   end
 end
 Agrad(:,jj)=Aa(:);
end

for ii=1:NElement
  ind=(Element(ii).Topology);               % The indices to g of the ii'th triangle.
  gg=g(ind,:);                              % A 3x2 or 6x2 matrix of triangle nodes in (x,y) coord.
  
 if any([Element(ii).Face{:,3}]),           %Checks if the triangle ii is the triangle that is
                                            % under the electrode.
    [In,Jn,InE]=find([Element(ii).Face{:,3}]);
    bind=Element(ii).Face{Jn,1};            % Nodes on the boundary
    ab=g(bind(:),:);

    if max(size(bind))==2                   % First order basis
     bb1=bound1([ab]);Bb1=zeros(max(size(ind)),1);
     bb2=bound2([ab]);Bb2=zeros(max(size(ind)));

     s(InE)=s(InE)+1/z(InE)*2*bb1; % 2*bb1 = length of the electrode.

     eind=[find(bind(1)==ind),find(bind(2)==ind)];
    else                                    % Second order basis
      bb1=boundquad1([ab]);Bb1=zeros(max(size(ind)),1);
      bb2=boundquad2([ab]);Bb2=zeros(max(size(ind)));

      s(InE)=s(InE)+1/z(InE)*electrlen([ab]);

      eind=[find(bind(1)==ind),find(bind(2)==ind),find(bind(3)==ind)];
    end

    Bb1(eind)=bb1;
    M(ind,InE)=M(ind,InE)-1/z(InE)*Bb1;

    Bb2(eind,eind)=bb2;
    Kb(ind,ind)=Kb(ind,ind)+1/z(InE)*Bb2;
  else                                      %The triangle isn't under the electrode.
  end
end  

S=sparse(diag(s));

[II1,C]=Current(Nel,NNode,'adj');
C=C(:,1:Nel-1);                             % For the voltage reference
C=sparse(C(:,1:Nel-1));                             

S=C'*S*C;
M=M*C;







