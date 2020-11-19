function [g2,H2,Ind2] = RefineMesh(g1,H1)

%RefineMesh Refines a coarse mesh by adding new nodes in the triangle edges
% Function [g2,H2,Ind2] = RefineMesh(g1,H1);
% Refines a coarse mesh by adding new nodes in the triangle edges. 
%
% INPUT
%
% g1 = node coordiantes of the coarse mesh
% H1 = connectivity of the coarse mesh
%
% OUTPUT
%
% g2 = node coordinates of the dense mesh
% H2 = connectivity of the dense mesh
% Ind2 = mapping between the two meshes


% Aku Seppänen 27.9.1999
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: aoseppan@venda.uku.fi


H2 = zeros(4*size(H1,1),size(H1,2));
gg = zeros(2*size(g1,1),size(g1,2));
gg(1:size(g1,1),:)=g1;

for i = 1:size(H1,1)
   i11=H1(i,1);
   i22=H1(i,2);
   i33=H1(i,3);
   g11 = g1(i11,:);
   g22 = g1(i22,:);
   g33 = g1(i33,:);
   g12 = (g11 + g22)/2;
   g13 = (g11 + g33)/2;
   g23 = (g22 + g33)/2;
   i12=size(g1,1)+(i-1)*3+1;
   i13=size(g1,1)+(i-1)*3+2;
   i23=size(g1,1)+i*3;
   gg([i12,i13,i23],:)=[g12;g13;g23];
   H2((i-1)*4+1:i*4,:)=[i11,i12,i13; i12,i22,i23; i13,i23,i33; i12,i23,i13];
end;

g2=gg(1,:);
for k=2:length(gg)
     I=find((g2(:,1)==gg(k,1)) & (g2(:,2)==gg(k,2)));
     J=find(H2(:)==k);
     if isempty(I)
        g2 = [g2; gg(k,:)];
        H2(J) = length(g2);
     else
        H2(J)=I;
     end;
end;

Ind2 = sparse(length(H2),length(H1));
for j=1:size(Ind2,2)
     Ind2((j-1)*4+1:j*4,j)=1;
end;



