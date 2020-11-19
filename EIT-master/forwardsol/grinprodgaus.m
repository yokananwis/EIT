function int=grinprodgaus(g,sigma);

%grinprodgaus Computes the integral of the product of the gradients in 2D FEM
% The function int=grinprodgaus(g,sigma);
% calculates the gradient part in the linear 
% FEM in 2D EIT.
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
L=[-1 1 0;-1 0 1];
Jt=L*g;
iJt=inv(Jt);
dJt=abs(det(Jt));
G=iJt*L;
%int=0;
% for ii=1:3
%  int=int+w(ii)*G'*G;
% end
%int=sigma*int*dJt;
int=1/2*sigma*G'*G*dJt;


