%EidorsDemoNodeSigmaquad Demonstrates the use of 2D EIT Package with real measurements, potential in quadratic basis
% EidorsDemoNodeSigmaquad Demonstrates the use of 2D EIT Package with real measurement
% using adjacent current injection patterns. The reconstructed image presents
% static (absolute) resistivity values.

% M. Vauhkonen 28.3.2000,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

load meshdata2nd  % Data for two different meshes.

NNode1=max(size(Node1));                      %The number of nodes
NElement1=max(size(Element1));                %The number of elements
NNode=max(size(Node));                      %The number of nodes
NElement=max(size(Element));                %The number of elements

H=reshape([Element.Topology],3,NElement)';
g=reshape([Node.Coordinate],2,NNode)';


L=16;                                     % The number of electrodes.
z=0.001*ones(L,1);                        % Contact impedances.
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
A=UpdateFemMatrix(Agrad,Kb,M,S,ones(size(Agrad,2),1));  % The system matrix.
[Uref,p,r]=ForwardSolution(NNode1,NElement1,A,C,T,MeasPatt,'real');
Urefel=Uref.Electrode(:);
Urefel=RemoveVolt(Urefel,L);

rho0=Urefel(:)\Uel(:);
sigma=1./rho0*ones(size(Agrad,2),1);


A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.
Uref=ForwardSolution(NNode1,NElement1,A,C,T,MeasPatt,'real',p,r);
Urefel=Uref.Electrode(:);
Urefel=RemoveVolt(Urefel,L);

J=JacobianNode(Node1,Element1,Agrad,Uref.Current,Uref.MeasField, ...
           [],'real');
J=RemoveJacob(J,L);
R=MakeRegmatrixNode(Node);
alpha=1.8;

W=diag(1./(Uel(:).^2));
iter=5;

no=[];
for ii=1:iter
 ii
 sigma=sigma+(J'*W*J+alpha*R'*R)\(J'*W*(Uel(:)-Urefel(:))-alpha*R'*R*sigma);
 no=[no;norm(Uel(:)-Urefel(:))^2+alpha*norm(R*sigma)^2];
 A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.
 Uref=ForwardSolution(NNode1,NElement1,A,C,T,MeasPatt,'real',p,r);
 Urefel=Uref.Electrode;
 Urefel=RemoveVolt(Urefel,L);
 J=JacobianNode(Node1,Element1,Agrad,Uref.Current,Uref.MeasField,[],'real');
 J=RemoveJacob(J,L);
 clf,Plotinvsolnode(sigma,g,H);drawnow;
end








