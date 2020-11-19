function [g,gp,H,EE]=cirgrid_eit(r,N,eI);

%cirgrid_eit Makes circular meshes for 2D EIT
% function [g,gp,H,EE]=cirgrid_eit(r,N,eI);
% g includes the grid points (nodes) and H the indices of 
% the triangles. See plcigrid to see the plotting and how to 
% get the true (rectangular) coordinates of the triangle points.
% r must be in descending order and must end with 0; the last 
% term in N must be 1. gp contains the polar coordinates of
% the nodes and E contains the indices of the elements that
% are on the boundary under the electrodes. eI is 2*1 vector,
% which contains the number of the elements under the electrode
% and number of the elements between the electrodes.
%
% INPUT
%
% r = vector of the radia
% N = vector of the nodes on each radius
% eI = element under and between the electrodes
%
% OUTPUT
%
% g = coordinates of the nodes (x,y)
% gp = polar coordinates of the nodes 
% H = element connectivity
% EE = elements under the electrodes having two nodes on the boundary

% J. Kaipio 12.4.1994. The innermost layer corrected 29.4.1994.
% Interior point indexes added 18.6.1994
% Modified by M. Vauhkonen 5.9.1994 from the version circgrid.m
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


nN=length(N);
e1=eI(1);e2=eI(2);
Nel=N(1)/(e1+e2); %The number of electrodes
phi=zeros(max(N)+1,nN);
g=[];
H=[];
E=[];
for ii=1:nN,
  phi(1:N(ii)+1,ii)=2*pi*[[0:N(ii)]]'/N(ii);
end
%phi(1:2,nN)=[0 0]';

for ii=1:nN
  for ij=1:N(ii)
    g=[g;[r(ii) phi(ij,ii)]];
  end
end
gp=g; % polar coordinates


% The first N(nN) points in g obviously not interior points
Int=[N(1)+1:sum(N)]';

for ii=1:nN-2
  icprev=1;
  for ij=1:N(ii)
      % Nearest point on inner layer
    cpt=abs(phi(1:N(ii+1)+1,ii+1) - (phi(ij,ii)+phi(ij+1,ii))/2);
    icpt=find(cpt==min(cpt));
    icpt=icpt(1);
      % Full circle termination settings
    if icpt==N(ii+1)+1,icpt1=1;else,icpt1=icpt;end
    if ij==N(ii),ij1=1;else,ij1=ij+1;end
    ipt1=find(g(:,1)==r(ii) & g(:,2)==phi(ij,ii));
    ipt2=find(g(:,1)==r(ii) & g(:,2)==phi(ij1,ii));
    ipt3=find(g(:,1)==r(ii+1) & g(:,2)==phi(icpt1,ii+1));
    H=[H;[ipt1 ipt2 ipt3]];
    % In E, there are the indices of boundary elements with two nodes
    % on the boundary 
    if ii==1,
      [Hr,Hc]=size(H);
      E=[E;Hr];
    end
% if ii==1, Hsize=size(H); Iborder=[Iborder;Hsize(1)];end
    for il=icprev:icpt-1
      ipt1=find(g(:,1)==r(ii) & g(:,2)==phi(ij,ii));
      ipt2=find(g(:,1)==r(ii+1) & g(:,2)==phi(il,ii+1));
      if il==N(ii+1),il1=1;else,il1=il+1;end
      ipt3=find(g(:,1)==r(ii+1) & g(:,2)==phi(il1,ii+1));
      H=[H;[ipt1 ipt2 ipt3]];
    end
      % Full circle termination
    if ij==N(ii)
      for il=icpt:N(ii+1)
          % Full circle termination settings
        if il==N(ii+1),il1=1;else il1=il+1;end
          %
        ipt1=find(g(:,1)==r(ii) & g(:,2)==phi(1,ii));
        ipt2=find(g(:,1)==r(ii+1) & g(:,2)==phi(il,ii+1));
        ipt3=find(g(:,1)==r(ii+1) & g(:,2)==phi(il1,ii+1));
        H=[H;[ipt1 ipt2 ipt3]];
      end
    end
    icprev=icpt;
  end
end

% The innermost triangles:
for ij=1:N(nN-1)-1
  ipt1=find(g(:,1)==r(nN-1) & g(:,2)==phi(ij,nN-1));
  ipt2=find(g(:,1)==r(nN-1) & g(:,2)==phi(ij+1,nN-1));
  ipt3=find(g(:,1)==0 & g(:,2)==0);
  H=[H;[ipt1 ipt2 ipt3]];
end
% Termination
ipt1=find(g(:,1)==r(nN-1) & g(:,2)==phi(N(nN-1),nN-1));
ipt2=find(g(:,1)==r(nN-1) & g(:,2)==phi(1,nN-1));
ipt3=find(g(:,1)==0 & g(:,2)==0);
H=[H;[ipt1 ipt2 ipt3]];

% Rectangular coordinates
g=[g(:,1).*cos(g(:,2)) g(:,1).*sin(g(:,2))];

% EE contains the elements under the electrodes

EE=zeros(Nel,e1);
NN=max(size(E));
ik=1;
  for ii=1:(e1+e2):NN-e1,
    EE(ik,:)=E(ii:(ii-1)+e1)';
    ik=ik+1;
  end























