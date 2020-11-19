function R=MakeRegmatrixFace(Element);

% Function R=MakeRegmatrixFace(Element);
% computes a regularization matrix R which is a
% version of 2D difference. 
%
% INPUT
%
% Element = element structure
%
% OUTPUT
%
% R = regularisation matrix

% M. Vauhkonen 13.8.1999, Univeristy of Kuopio, Finland

Facelist=[];
for ii=1:max(size(Element))
 Fl=[Element(ii).Face{:,2}]';
 fFl=find(Fl);
 add=[ones(max(size(fFl)),1)*ii,Fl(fFl)];
 Facelist=[Facelist;add];
end
Facelist=sort(Facelist')';
Facelist=unique(Facelist,'rows');
[rN,cN]=size(Facelist);
R=sparse(rN,max(max(Facelist)));
for ii=1:rN
 R(ii,Facelist(ii,:))=[1,-1];
end
  
