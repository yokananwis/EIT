%EidorsDemoNodeSigma Demonstrates the use of 2D EIT Package with real measurements, potential in linear basis
% EidorsDemoNodeSigma Demonstrates the use of 2D EIT Package with real measurement
% using adjacent current injection patterns. The reconstructed image presents
% static (absolute) conductivity values.

% M. Vauhkonen 28.3.2000,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

load meshdata  % Data for two different meshes.

NNode1=max(size(Node1));                      %The number of nodes
NElement1=max(size(Element1));                %The number of elements
g=reshape([Node1.Coordinate],2,NNode1)';
H=reshape([Element1.Topology],3,NElement1)';



L=16;                                     % The number of electrodes.
z=0.005*ones(L,1);                        % Contact impedances.
[II1,T]=Current(L,NNode1,'adj');          % Adjacent current pattern.

%Load some data
load bubble2.dat
meas=bubble2(1280-255:1280);
Uel=reshape(meas,16,16);
Uel=RemoveVolt(Uel,L);

%%             PROCEDURE TO SOLVE THE INVERSE PROBLEM           %%

% Approximate the best homogenous resistivity.

MeasPatt=-T;
[Agrad,Kb,M,S,C]=FemMatrixNode(Node1,Element1,z);
A=UpdateFemMatrix(Agrad,Kb,M,S,ones(NNode1,1));  % The system matrix.
[Uref,p,r]=ForwardSolution(NNode1,NElement1,A,C,T,MeasPatt,'real');
Urefel=Uref.Electrode(:);
Urefel=RemoveVolt(Urefel,L);

rho0=Urefel(:)\Uel(:);
sigma=1./(rho0*ones(NNode1,1));


A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.
Uref=ForwardSolution(NNode1,NElement1,A,C,T,MeasPatt,'real',p,r);
Urefel=Uref.Electrode(:);
Urefel=RemoveVolt(Urefel,L);

J=JacobianNode(Node1,Element1,Agrad,Uref.Current,Uref.MeasField, ...
           [],'real');
J=RemoveJacob(J,L);
R=MakeRegmatrixNode(Node1);
alpha=0.5;

iter=10;

for ii=1:iter
 sigma=sigma+(J'*J+alpha*R'*R)\(J'*(Uel(:)-Urefel(:))-alpha*R'*R*sigma);
 A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.
 Uref=ForwardSolution(NNode1,NElement1,A,C,T,MeasPatt,'real',p,r);
 Urefel=Uref.Electrode;
 Urefel=RemoveVolt(Urefel,L);
 J=JacobianNode(Node1,Element1,Agrad,Uref.Current,Uref.MeasField,[],'real');
 J=RemoveJacob(J,L);
 clf,Plotinvsolnode(sigma,g,H);drawnow;
end








