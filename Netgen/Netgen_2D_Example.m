cyl_shp = [0, 10, 0.5];
[fmdl,mat_idx]= ng_mk_cyl_models(cyl_shp,16,0.02);
figure(1);
title("Finite Element Model");
show_fem(fmdl);
% show_fem(fmdl, [0 1.015 0]);

% Nod=fmdl.nodes;
% Topol=fmdl.elems;
% clf,plcigrid(Nod,Topol,[],'r'),axis equal;
% NTopol=max(size(Topol));
% x = [4 8];
% y = [4 8];
% r=norm([x(1),y(1)]-[x(2),y(2)]);
% Ind = sparse(zeros(NTopol,1));
% 
% for ij=1:NTopol
%  nodes=Nod(Topol(ij,:),:);
%  center = mean(nodes);
%     if norm(center(:)-[x(1);y(1)])<r
%        Ind(ij) = 1;
%     end
% end
% Ind=find(Ind);
%  for ii=1:max(size(Ind))
%   Hii=Nod(Topol(Ind(ii),:),:);
%   patch(Hii(:,1),Hii(:,2),1);
%  end
% 
% Rh=675;
% Rn=17500;
% disp('Generating Finite Element Model.')
% Ind=Set_Resistivity(fmdl.nodes,fmdl.elems);       % Make data for an inhomogeneity.
% sigma=1/Rh*ones(max(size(fmdl.elems)),1);            % Make a conductivity vector.
% sigma(Ind)=1/Rn;			  % Conductivity of the inhomogeneity.
% % sigma = CreateInhomogeneities(Node2,Element2,7);
% 
% % Eventually we'll want to get rid of Plotinvsol or rewrite it.*
% figure(2)
% clf,Plotinvsol(1./sigma,fmdl.nodes,fmdl.elems);colorbar,title('Finite Element Model');