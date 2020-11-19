function int=electrlen(g);

%electrlen Computes the length of the electrode in the case of isoparametric (quadratic) triangular element 
% Function int=electrlen(g) computes the length of the electrode in the case of isoparametric triangular element
% The electrode is the curve defined by the coordinates in g.
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


w=[1/2,1/2];
ip=[1/2-1/6*sqrt(3),1/2+1/6*sqrt(3)];
int=0;
 for ii=1:2
  dJt=sqrt((g(1,1)*(4*ip(ii)-3)+g(2,1)*(4-8*ip(ii))+g(3,1)*(4*ip(ii)-1))^2+ ...
            (g(1,2)*(4*ip(ii)-3)+g(2,2)*(4-8*ip(ii))+g(3,2)*(4*ip(ii)-1))^2);
  int=int+w(ii)*dJt;
 end










