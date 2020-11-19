function int=bound1(g);

%bound1 Computes the boundary integral of one linear basis function in 2D FEM
% Function int=bound1(g) calculates the boundary integral
% of one basis function from g(1,:) to g(2,:).
%
% INPUT
%
% g = integration endpoints
%
% OUTPUT
%
% int = value of the integral

% 10.5. 1996 P. Ronkanen and M. Vauhkonen
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


int=1/2*sqrt((g(2,1)-g(1,1))^2+(g(2,2)-g(1,2))^2); 

%w=[1/2,1/2];
%ip=[1/2-1/6*sqrt(3),1/2+1/6*sqrt(3)];
%dJt=sqrt((g(2,1)-g(1,1))^2+(g(2,2)-g(1,2))^2); 
%int=0;
% for ii=1:2
%  S=[1-ip(ii);ip(ii)];
%  int=int+w(ii)*S;
% end
%int=int*dJt;

