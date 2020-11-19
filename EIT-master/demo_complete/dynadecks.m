function []=dynadecks(X,WW,g,H,nc);

%dynadecks Plots the results of the dynamic reconstruction computed by Dyneit.m  
% function []=dynadecks(X,WW,g,H,nc);
% Plots the results of the dynamic reconstruction computed by Dyneit.m
%
% INPUT
%
% X = matrix of the reconstructions
% WW = mapping from sparse to dense mesh (optional, use []) 
% g = coordinates of the mesh nodes
% H = connectivity data
% nc = the number of columns in the produced figure

% Pasi A. Karjalainen, 1997
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Pasi.Karjalainen@uku.fi

[N,n]=size(X);
N=sqrt(N);
mx=max(max(X));
minx=min(min(X));

if nargin<2,
 nc=floor(sqrt(n));
end

nr=nc+ceil((n/nc)-nc);
eps=0.00001;
voffs=0.005;
hoffs=0.005;
colsep=(1/(nc)); 
rowsep=(1/(nr)); 
wdt=colsep-1*hoffs;
hgt=rowsep-1*voffs;
i=1:n;
row=nr-1-floor(i/nc-1/nc+eps);
col=fix(nc*((i/nc-1/nc+eps)-floor(i/nc-1/nc+eps)));
for i=1:n,
ax(i)=axes;
set(ax(i),'position',[col(i)*colsep+hoffs row(i)*rowsep+voffs wdt hgt]);

if ~isempty(WW)
 Plotinvsol(WW*X(:,i),g,H),set(gca,'climmode','manual')
 set(gca,'clim',[minx,mx]),axis('image')
else
 Plotinvsol(X(:,i),g,H),set(gca,'climmode','manual')
 set(gca,'clim',[minx,mx]),axis('image')
end
end




