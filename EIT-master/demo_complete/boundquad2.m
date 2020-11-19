function int=boundquad2(g);

%boundquad2 Computes the boundary integral of the product of two quadratic basis functions 
% Function int=boundquad2(g) calculates the boundary integral
% of the product of two quadratic basis function over the curve defined by the coordinates in g.
%
% INPUT
%
% g = global coordinates of the integration curve
%
% OUTPUT
%
% int = value of the integral

% 10.5. 1996 P. Ronkanen and M. Vauhkonen
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

w=[5/18,8/18,5/18];
ip=[1/2-1/10*sqrt(15),1/2,1/2+1/10*sqrt(15)];
int=0;
 for ii=1:3
  S=[2*ip(ii)^2-3*ip(ii)+1; ...
     -4*ip(ii)^2+4*ip(ii); ...
     2*ip(ii)^2-ip(ii)];
  dJt=sqrt((g(1,1)*(4*ip(ii)-3)+g(2,1)*(4-8*ip(ii))+g(3,1)*(4*ip(ii)-1))^2+ ...
            (g(1,2)*(4*ip(ii)-3)+g(2,2)*(4-8*ip(ii))+g(3,2)*(4*ip(ii)-1))^2);
  int=int+w(ii)*S*S'*dJt;
 end















