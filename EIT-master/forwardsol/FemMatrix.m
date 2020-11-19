function [Agrad,Kb,M,S,C]=FemMatrix(Node,Element,z);

%FemMatrix Computes the blocks of the system matrix for 2D EIT with linear and quadratic basis
% Function [Agrad,Kb,M,S,C]=FemMatrix(Node,Element,z); 
% computes the matrices needed in the finite element 
% approximation of the 2D EIT forward problem. 
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

% M. Vauhkonen 11.5.1994, modified from the version of J.P. Kaipio
% 25.4.1994. Modified 5.9.1994 by M. Vauhkonen for EIT.
% Modified 13.8.1999 and 23.3.2000 for the EIDORS by M. Vauhkonen, 
% University of Kuopio, Department of Applied Physics, PO Box 1627, 
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi 

Nel=max(size(z));                           %The number of electrodes.
NNode=max(size(Node));                      %The number of nodes
NElement=max(size(Element));                %The number of elements
M=sparse(NNode,Nel);
Kb=sparse(NNode,NNode);
Agrad=sparse(NNode^2,NElement);
s=zeros(Nel,1);
g=reshape([Node.Coordinate],2,NNode)';      %Nodes

for ii=1:NElement 
  A=sparse(NNode,NNode);
  ind=(Element(ii).Topology);               % Indices f the ii'th element
  gg=g(ind,:);                              % A 3x2 or 6x2 matrix of triangle vertices in (x,y) coord.
  if max(size(gg))==3                       % First order basis   
   grint=grinprodgaus(gg,1);
  else
   grint=grinprodgausquad(gg,1);            % Second order basis
  end

 if any([Element(ii).Face{:,3}]),           % Checks if the triangle ii is a triangle that is
                                            % under the electrode.
    [In,Jn,InE]=find([Element(ii).Face{:,3}]); 
    bind=Element(ii).Face{Jn,1};            % Nodes on the boundary
    ab=g(bind(:),:);

    if max(size(bind))==2                   % First order basis
     bb1=bound1([ab]);Bb1=zeros(max(size(ind)),1);
     bb2=bound2([ab]);Bb2=zeros(max(size(ind)));

     s(InE)=s(InE)+1/z(InE)*2*bb1; % 2*bb1 = length of the electrode.

     eind=[find(bind(1)==ind),find(bind(2)==ind)];   
    else				    % Second order basis
      bb1=boundquad1([ab]);Bb1=zeros(max(size(ind)),1);
      bb2=boundquad2([ab]);Bb2=zeros(max(size(ind)));

      s(InE)=s(InE)+1/z(InE)*electrlen([ab]); 

      eind=[find(bind(1)==ind),find(bind(2)==ind),find(bind(3)==ind)];
    end 

    Bb1(eind)=bb1;   
    M(ind,InE)=M(ind,InE)-1/z(InE)*Bb1;
    
    Bb2(eind,eind)=bb2;
    A(ind,ind)=grint;
    Agrad(:,ii)=A(:);
    Kb(ind,ind)=Kb(ind,ind)+1/z(InE)*Bb2;    
 else                                       % The triangle isn't under the electrode.
    A(ind,ind) = grint; 
    Agrad(:,ii)=A(:); 
 end
end  

S=sparse(diag(s));


[II1,C]=Current(Nel,NNode,'adj');
C=C(:,1:Nel-1);                             % For the voltage reference
C=sparse(C(:,1:Nel-1));                             

S=C'*S*C;
M=M*C;







