function [I,T]=Current(L,lg,style,rms,numpat);

%Current Supplies some current patterns for 2D EIT
% Function [I,T]=Current(L,lg,style,rms,numpat) calculates 
% different current patterns. The 'style' can be 'tri' (trigonometric), 
% 'adj' (adjacent), 'opp' (opposite) or 'ref' (one reference). 
%
% INPUT
%
% L = the number of electrodes
% lg = the number of nodes
% style = 'tri' (trigonometric), 'adj' (adjacent), 'opp' (opposite) or 'ref' (one reference)
% rms = the root mean square value of the injected current. Default=1 (optional).
% numpat = the number of patterns (optional, in certain cases forced to be L-1)
%
% OUTPUT
%
% I = the total current matrix, internal and on the electrodes 
% T = currents on the electrodes

% M. Vauhkonen 6.9.1994. Modified 13.9.1994 for NOSER.
% Modified 13.8.1999 by M. Vauhkonen,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


if nargin < 4 & style == 'tri'
 rms=1;
 numpat=L-1;
end

if nargin < 4 & style ~= 'tri'
 rms=1;
 numpat=L;
end

if nargin == 4 & style == 'tri'
 numpat=L-1;
end

if nargin == 4 & style ~= 'tri'
 numpat=L;
end


switch style
case 'tri'
if numpat > L-1,
  str=['Too many current patterns. Automatically set to ', num2str(L-1),'.' ];
  disp(str)
  numpat =L-1;
 end
II=zeros(L,numpat);
l=(1:L)';
th=2*pi*l/L;
 for k=1:(numpat+1)/2,
  II(:,k)=rms*cos(k*th);
 end
 for k=(numpat+1)/2+1:numpat
  II(:,k)=rms*sin((k-L/2)*th);
 end
I=[zeros(lg,numpat);II];
T=II;

case 'adj'
II=zeros(L,numpat);
l=(1:L)';
II1=diag(ones(L,1));
II2=diag(ones(L-1,1),-1);
if numpat < L
 II=II1(:,1:numpat)-II2(:,1:numpat);
else
 II=II1(:,1:numpat-1)-II2(:,1:numpat-1);
 II=[II,[-1;zeros(L-2,1);1]];% Optional!!
end
I=[zeros(lg,numpat);rms*II];
T=II;

case 'ref'
 if numpat > L-1, 
  str=['Too many current patterns. Automatically set to ', num2str(L-1),'.' ];
  disp(str)
  numpat =L-1;
 end
II=zeros(L,numpat);
II1=diag(ones(L-1,1),-1);
II(1,:)=ones(1,numpat);
T=II-II1(:,1:numpat);
I=[zeros(lg,numpat);rms*T];


case 'opp'
if numpat > L/2, 
str=['Too many current patterns. Automatically set to ', num2str(L/2),'.' ];
  disp(str)
numpat=L/2;end
II=zeros(L,numpat);
l=(1:L)';
II1=diag(ones(L,1));
II2=diag(ones(L,1),-L/2);
II=II1(:,1:numpat)-II2(1:L,1:numpat);
%II=[II(:,1:ceil(numpat/2)),-II(:,1:floor(numpat/2))];% Optional!!
I=[zeros(lg,numpat);rms*II];
T=II;
end











