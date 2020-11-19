function [U,p,r]=ForwardSolution(NNode,NElement,A,C,T,MeasPattern,style,p,r)

%ForwardSolution Solves the potential distribution in 2D EIT 
% Function [U,p,r]=ForwardSolution(NNode,NElement,A,C,T,MeasPattern,style,p,r);
% computes the internal voltages for the injected currents (U.Current)
% and for the "measurement field" (U.MeasField) and the voltages
% on the electrodes (U.Electrode).
%
% INPUT
%
% NNode = the number of nodes
% NElement = the number of elements
% A = the system matrix
% C = the reference matrix
% T = current patterns
% MeasPattern = measurement pattern, default \sum_{l=1}^L=0.
% style = 'comp' for complex reconstruction and 'real' for real
% p = the permutation vector (optional)
% r = reversed permutation   (optional)
%
% OUTPUT
%
% U = potential data structure, includes potential corresponding to the actual
%     current injections (U.Current), for the measurement pattern (U.MeasField)
%     and voltages on the electrodes (U.Electrode)
% p = the permutation vector (optional)
% r = reversed permutation   (optional)
         
% M. Vauhkonen 13.8.1999
% University of Kuopio, Department of Applied Physics,
% PO Box 1627, FIN-70211, Kuopio, Finland, email:Marko.Vauhkonen@uku.fi
% Ordering added 15.11.1999, modified 28.3.2000 by M. Vauhkonen

L=max(size(A))-NNode+1;    % The number of electrodes     

if strcmp(style,'comp')
 II1=sparse([zeros(L,NNode),C,zeros(L,NNode+L-1);zeros(L,2*NNode+L-1),C]);
  if ~isempty(MeasPattern) 
   II1=MeasPattern'*II1;
   II1=II1'; 
  else
   MeasPattern=eye(2*max(size(C)));
   II1=II1';
  end
 II=sparse([[zeros(NNode,size(T,2));C'*T];zeros(NNode+L-1,size(T,2))]);
 A=[real(A),-imag(A);imag(A),real(A)];
 UU=A'\II1; %Voltages for the "measurement field" and for the current patterns.
 UU=full([UU,A\II]);
 U.MeasField=[UU(1:NNode,1:size(II1,2));UU(NNode+L:2*NNode+L-1,1:size(II1,2))]; %The "measurement field" data
 U.Electrode=MeasPattern'*[C,zeros(size(C));zeros(size(C)),C]* ...
            [UU(NNode+1:NNode+L-1,size(II1,2)+1:size(UU,2));UU(2*NNode+L:size(UU,1), ...
            size(II1,2)+1:size(UU,2))]; %Voltages on the electrodes
 U.Current=[UU(1:NNode,size(II1,2)+1:size(UU,2));UU(NNode+L:2*NNode+L-1,size(II1,2)+1:size(UU,2))];
 
 if nargin<8
  p=symamd(A);
  r(p)=1:max(size(p));
 end

elseif strcmp(style,'real')
 II1=sparse([zeros(L,NNode),C]);
  if ~isempty(MeasPattern) 
   II1=MeasPattern'*II1;
   II1=II1'; 
  else
   MeasPattern=eye(max(size(C)));
   II1=II1';
  end
 II=[zeros(NNode,size(T,2));C'*T];

 if nargin<8
  p=symamd(A);
  r(p)=1:max(size(p));
 end

 % I'm honestly not sure what the hell is happening here.
 R=chol(A(p,p));
 UU=R\(R'\[II1(p,:),II(p,:)]);
 UU=full(UU(r,:));
 U.MeasField=UU(1:NNode,1:size(II1,2)); %The "measurement field" data
 U.Electrode=MeasPattern'*C*UU(NNode+1:size(A,1),size(II1,2)+1:size(UU,2));%Voltages on the electrodes
 U.Current=UU(1:NNode,size(II1,2)+1:size(UU,2));
end















