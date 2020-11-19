%EidorsDemo2quad Demonstrates the use of 2D EIT Package with real measurements, potential in quadratic basis
% EidorsDemo2quad Demonstrates the use of 2D EIT Package with real measurements  
% using adjacent current injection patterns. The reconstructed image presents 
% static (absolute) resistivity values. 
% See also EidorsDemo2

% M. Vauhkonen 28.3.2000,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

load meshdata2nd  % Data for two different meshes.

NNode1=max(size(Node));                      %The number of nodes
NElement1=max(size(Element));                %The number of element
NNode2=max(size(Node1));                      %The number of nodes
NElement2=max(size(Element1));                %The number of elements

g1=reshape([Node.Coordinate],2,NNode1)';
H1=reshape([Element.Topology],3,NElement1)';
g2=reshape([Node1.Coordinate],2,NNode2)';
H2=reshape([Element1.Topology],6,NElement2)';


L=16;					  % The number of electrodes.
z=0.005*ones(L,1);			  % Contact impedances.
[II1,T]=Current(L,NNode2,'adj');	  % Adjacent current pattern.

%Load some data
load bubble2.dat
meas=bubble2(1280-255:1280);
Uel=reshape(meas,16,16);
Uel=RemoveVolt(Uel,L);

%%             PROCEDURE TO SOLVE THE INVERSE PROBLEM           %%

% Approximate the best homogenous resistivity.

MeasPatt=-T;
[Agrad,Kb,M,S,C]=FemMatrix(Node1,Element1,z);
A=UpdateFemMatrix(Agrad,Kb,M,S,ones(NElement2,1));  % The system matrix.
[Uref,p,r]=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real');
Urefel=Uref.Electrode;
Urefel=RemoveVolt(Urefel,L);

rho0=Urefel(:)\Uel(:);
rho=rho0*ones(NElement2,1);

A=UpdateFemMatrix(Agrad,Kb,M,S,1./rho);  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real',p,r);
Urefel=Uref.Electrode(:);
Urefel=RemoveVolt(Urefel,L);

J=Jacobian(Node1,Element1,Agrad,Uref.Current,Uref.MeasField, ...
           rho,'real');
J=RemoveJacob(J,L);

% Regularisation parameter and matrix
alpha=0.001; 
R=MakeRegmatrix(Element1);

iter=6;

for ii=1:iter
 rho=rho+(J'*J+alpha*R'*R)\(J'*(Uel(:)-Urefel(:))-alpha*R'*R*rho);
 A=UpdateFemMatrix(Agrad,Kb,M,S,1./rho);  % The system matrix.
 Uref=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real',p,r);
 Urefel=Uref.Electrode;
 Urefel=RemoveVolt(Urefel,L);
 J=Jacobian(Node1,Element1,Agrad,Uref.Current,Uref.MeasField,rho,'real');
 J=RemoveJacob(J,L);
 clf,Plotinvsol(rho,g1,H1);drawnow;
end






