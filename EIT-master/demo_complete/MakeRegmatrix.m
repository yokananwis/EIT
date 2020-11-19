function R=MakeRegmatrix(Element);

%MakeRegmatrix Computes a regularisation matrix which includes smoothness assumptions
% Function R=MakeRegmatrix(Element);
% computes a regularization matrix R which is a
% version of 2D difference. 
%
% INPUT
% Element = element structure
%
% OUTPUT
%
% R = regularisation matrix

% M. Vauhkonen 13.8.1999, 
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


msE=max(size(Element));
R=sparse(msE,msE);
for ii=1:msE
 Inds=[Element(ii).Face{:,2}];
 R(ii,Inds(Inds~=0))=-1;
 R(ii,ii)=abs(sum(R(ii,:)));
end
  


