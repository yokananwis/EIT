function Vh=RemoveVolt(Uel,L)

% I DON'T THINK THIS CURRENTLY WORKS IF NUMPAT < L-1.

%RemoveVolt Removes all the voltages measured on the current carrying electrodes
% Function Vh=RemoveVolt(Uel,L);
% removes all the voltages measured on the
% current carrying electrodes.
% Change the Jacobian accordingly, see RemoveJacob.m.
%
% INPUT
%
% Uel = voltages on the electrodes
% L = number of electrodes
% 
% OUTPUT
%
% Vh = new voltages 


% M. Vauhkonen 6.9.1999, 
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

  Uel=Uel(:);
  [rUel,~]=size(Uel);
  Uh=reshape(Uel,L,rUel/L); % Might need some rounding here.
  Vadj=Uh;
  kk=1;
   for ii=1:L
    Vadj(:,ii)=[Vadj(kk:L,ii);Vadj(1:kk-1,ii)];
    kk=kk+1;
   end
  Vadj([1 2 L],:)=[];
  Vh=Vadj;
  



