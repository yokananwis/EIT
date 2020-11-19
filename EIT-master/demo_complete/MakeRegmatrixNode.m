function R=MakeRegmatrixNode(Node);

%MakeRegmatrixNode Computes a regularisation matrix which includes smoothness assumptions. Conductivity in linear basis.
% Function R=MakeRegmatrixNode(Node);
% computes a regularization matrix R which is an
% version of 2D difference. Conductivity is approximated
% with linear basis.
%
% INPUT
%
% Node = nodal structure
%
% OUTPUT
%
% R = regularisation matrix
%
% M. Vauhkonen 13.8.1999,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


msN=max(size(Node));
R=sparse(msN,msN);
for ii=1:msN
 Inds=[Node(ii).NodeConnection];
 R(ii,Inds(Inds~=ii))=-1;
 R(ii,ii)=abs(sum(R(ii,:)));
end
  


