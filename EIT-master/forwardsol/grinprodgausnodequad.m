function int=grinprodgausnodequad(g,I);

%grinprodgausnodequad Computes the gradient part of the system matrix
%in the quadratic basis case. The conductivity  is in the linear basis.
% Function int=grinprodgausnodequad(g,I);
% computes the gradient part of the system matrix
% in the quadratic basis case. The conductivity  is in the linear basis.
%
% INPUT
%
% g = nodes of the element
% I = index of  the chosen node 
%
% OUTPUT
%
% int = three dimensional array or 3x3 matrix of the values of integral

% P. Ronkanen and M. Vauhkonen 10.5. 1996. Modified by M. Vauhkonen for EIDORS
% 11.5.2000.
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

w=[1/6*ones(3,1)];
ip=[1/2 0;1/2 1/2;0 1/2];

int=0;
 for ii=1:3
  S=[1-ip(ii,1)-ip(ii,2);ip(ii,1);ip(ii,2)];
  L=[4*(ip(ii,1)+ip(ii,2))-3, -8*ip(ii,1)-4*ip(ii,2)+4, ...
     4*ip(ii,1)-1, 4*ip(ii,2), 0, -4*ip(ii,2); ...
     4*(ip(ii,1)+ip(ii,2))-3, -4*ip(ii,1), ...
     0, 4*ip(ii,1), 4*ip(ii,2)-1, -8*ip(ii,2)-4*ip(ii,1)+4];
  Jt=L*g;
  iJt=inv(Jt);
  dJt=abs(det(Jt));
  G=iJt*L;
  int=int+w(ii)*S(I)*G'*G*dJt;
 end







