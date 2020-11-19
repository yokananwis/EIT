function int=grinprodgausquad(g,sigma);

%grinprodgausquad Computes the integral of the product of the gradients in 2D FEM for quadratic isoparametric triangular elements
% Function int=grinprodgausquad(g,sigma);
% calculates the integral of the product of the gradients in 2D FEM for quadratic
% isoparametric triangular elements
%
%
% INPUT
%
% g = nodal coordinates
% sigma = conductivity of the element
%
% OUTPUT
%
% int = value of the integral

% P. Ronkanen and M. Vauhkonen 10.5. 1996
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


w=[1/6*ones(3,1)];
ip=[1/2 0;1/2 1/2;0 1/2];

int=0;
 for ii=1:3
 S=[2*ip(ii,1)^2+2*ip(ii,2)^2+4*ip(ii,1)*ip(ii,2)-3*ip(ii,1)-3*ip(ii,2)+1; ...
  -4*ip(ii,1)^2-4*ip(ii,1)*ip(ii,2)+4*ip(ii,1); ...
   2*ip(ii,1)^2-ip(ii,1);4*ip(ii,1)*ip(ii,2);2*ip(ii,2)^2-ip(ii,2); ...
  -4*ip(ii,2)^2-4*ip(ii,1)*ip(ii,2)+4*ip(ii,2)];
  L=[4*(ip(ii,1)+ip(ii,2))-3, -8*ip(ii,1)-4*ip(ii,2)+4, ...
     4*ip(ii,1)-1, 4*ip(ii,2), 0, -4*ip(ii,2); ...
     4*(ip(ii,1)+ip(ii,2))-3, -4*ip(ii,1), ...
     0, 4*ip(ii,1), 4*ip(ii,2)-1, -8*ip(ii,2)-4*ip(ii,1)+4];
  Jt=L*g;
  iJt=inv(Jt);
  dJt=abs(det(Jt));
  G=iJt*L;
  int=int+w(ii)*G'*G*dJt;
 end
int=sigma*int;



