function int=grinprodgausnode(g,I);

%grinprodgausnode Computes the gradient part of the system matrix 
%in the linear basis case. The conductivity  is also in the linear basis. 
% Function int=grinprodgausnode(g,I);
% computes the gradient part of the system matrix 
% in the linear basis case. The conductivity  is also in the linear basis. 
%
% INPUT
%
% g = nodes of the element
% I = index of  the chosen node (optional)
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
L=[-1 1 0;-1 0 1];
Jt=L*g;
iJt=inv(Jt);
dJt=abs(det(Jt));
G=iJt*L;
GdJt=G'*G*dJt;
msg=max(size(g));

if nargin > 1
 S=[[1-ip(1,1)-ip(1,2);ip(1,1);ip(1,2)],[1-ip(2,1)-ip(2,2);ip(2,1);ip(2,2)],[1-ip(3,1)-ip(3,2);ip(3,1);ip(3,2)]]*w;
 int=S(I)*GdJt;
else
  int=zeros(msg,msg,3);
  S=[[1-ip(1,1)-ip(1,2);ip(1,1);ip(1,2)],[1-ip(2,1)-ip(2,2);ip(2,1);ip(2,2)],[1-ip(3,1)-ip(3,2);ip(3,1);ip(3,2)]]*w;
 for jj=1:3
  int(:,:,jj)=S(jj)*GdJt;
 end
end






