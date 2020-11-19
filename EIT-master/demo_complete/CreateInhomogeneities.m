function [sigma] = CreateInhomogeneities(Node2,Element2,baseR)
% CREATEINHOMOGENEITIES will create a domain to replicate the cross section
% of lungs and heart. Returns a conductivity grid sigma. Must input the
% base resistivity, baseR.

% Assign resistivity values, measured in ohm-meters.
% Bone: 150
% Fat: 25
% Muscle: 5
% Blood: 1.5
% Plasma/interstitial fluid: 0.65
% The above values come from the upper arm.
% "Lung resistivity is about five-fold higher than other organs in the
% thorax."

% For now, choose background resistivity to be 8 ohm*m. Heart is full of
% blood and can be 1.5 ohm*m. Lungs can be 8*5=40 ohm*m.

%lungR = 5*baseR;
lungR = 35;
heartR = 1.5;

NElement2=max(size(Element2));                %The number of elements

lungInd1 = ChooseCircleXY(Node2,Element2,[-5,3.8],[-5,9]);
lungInd2 = ChooseCircleXY(Node2,Element2,[5,3.8],[5,9]);
heartInd = ChooseCircleXY(Node2,Element2,[0,-4],[0,-6]);
resis = baseR*ones(NElement2,1);         % Resistivity vector
resis(lungInd1)=lungR;
resis(lungInd2)=lungR;
resis(heartInd)=heartR;
sigma=1./resis;                        % Make a conductivity vector.

