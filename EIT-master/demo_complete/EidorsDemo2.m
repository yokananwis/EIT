% EidorsDemo2 Demonstrates the use of 2D EIT Package with real measurement 
% using adjacent current injection patterns. The reconstructed image presents 
% static (absolute) resistivity values. 

% M. Vauhkonen 28.3.2000,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

load meshdata  % Data for two different meshes.

NNode1=max(size(Node1));                      %The number of nodes
NElement1=max(size(Element1));                %The number of element
NNode2=max(size(Node2));                      %The number of nodes
NElement2=max(size(Element2));                %The number of elements

g1=reshape([Node1.Coordinate],2,NNode1)';
H1=reshape([Element1.Topology],3,NElement1)';
g2=reshape([Node2.Coordinate],2,NNode2)';
H2=reshape([Element2.Topology],3,NElement2)';


L=16;					  % The number of electrodes.
z=0.005*ones(L,1);			  % Contact impedances.
[II1,T]=Current(L,NNode2,'adj');	  % Adjacent current pattern.
% I don't know why we need adjacent pattern for this.

% Load some data measured with the EIT system built in Kuopio
% So we have to create this data ahead of time, maybe.
load bubble2.dat
meas=bubble2(1280-255:1280);
Uel=reshape(meas,L,L);
Uel=RemoveVolt(Uel,L);

%%             PROCEDURE TO SOLVE THE INVERSE PROBLEM           %%

% Approximate the best homogenous resistivity.

MeasPatt=-T;
[Agrad,Kb,M,S,C]=FemMatrix(Node2,Element2,z);
Agrad1=Agrad*Ind2;
A=UpdateFemMatrix(Agrad,Kb,M,S,ones(NElement2,1));  % The system matrix.
[Uref,p,r]=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real');
Urefel=Uref.Electrode;
Urefel=RemoveVolt(Urefel,L);

rho0=Urefel(:)\Uel(:);


A=UpdateFemMatrix(Agrad,Kb,M,S,1/rho0*ones(NElement2,1));  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real',p,r);
Urefel=Uref.Electrode(:);
Urefel=RemoveVolt(Urefel,L);

rho=rho0*ones(size(Ind2,2),1);
J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField, ...
           rho,'real');
J=RemoveJacob(J,L);

% Regularisation parameter and matrix

alpha=0.001; 
R=MakeRegmatrix(Element1);

iter=6;

no=norm(Uel(:)-Urefel(:))^2+alpha*norm(R*rho)^2;
for ii=1:iter
 rho=rho+(J'*J+alpha*(R')*R)\(J'*(Uel(:)-Urefel(:))-alpha*(R')*R*rho);
 no=[no;norm(Uel(:)-Urefel(:))^2+alpha*norm(R*rho)^2];             %Error norm
 rhobig=Ind2*rho;
 A=UpdateFemMatrix(Agrad,Kb,M,S,1./rhobig);  % The system matrix.
 Uref=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real',p,r);
 Urefel=Uref.Electrode;
 Urefel=RemoveVolt(Urefel,L);
 J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,rho,'real');
 J=RemoveJacob(J,L);
 %figure(ii)
 clf,Plotinvsol(rho,g1,H1);colorbar;drawnow;
end






