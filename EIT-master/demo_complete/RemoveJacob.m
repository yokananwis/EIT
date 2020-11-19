function Jh=RemoveJacob(J,L);

%RemoveJacob Removes the derivatives corresponding to the current  carrying electrodes
% Function Jh=RemoveJacob(J,L);
% removes all the derivatives that correspond to the
% current carrying electrodes.
% Change the voltages accordingly, see RemoveVolt.m.
%
% INPUT
%
% J = Jacobian
% L = number of electrodes
%
% OUTPUT
%
% Jh = new Jacobian

% M. Vauhkonen 6.9.1999, 
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi



[rJ,cJ]=size(J);
J=reshape(J,L,rJ/L*cJ);
kk=1;p=1;
for jj=1:cJ
 for ii=1:L
  J(:,p)=[J(kk:L,p);J(1:kk-1,p)];
 kk=kk+1;p=p+1;
end
kk=1;
end
J([1 2 L],:)=[];
Jh=reshape(J,(L-3)*L,cJ);

  

