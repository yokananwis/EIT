%EidorsDemoDynamic Demonstrates the use of 2D EIT Package for dynamic imaging
% EidorsDemoDynamic Demonstrates the use of 2D EIT Package for dynamic imaging.
% Here only a static image is reconstructed that shows the convergence of the algorithms.

% M. Vauhkonen 28.3.2000
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

load meshdata % Data for two different meshes.

NNode1=max(size(Node1));                      %The number of nodes
NElement1=max(size(Element1));                %The number of element
NNode2=max(size(Node2));                      %The number of nodes
NElement2=max(size(Element2));                %The number of elements

g1=reshape([Node1.Coordinate],2,NNode1)';
H1=reshape([Element1.Topology],3,NElement1)';
g2=reshape([Node2.Coordinate],2,NNode2)';
H2=reshape([Element2.Topology],3,NElement2)';

Ind=ChooseCircle(Node2,Element2);         % Make data for a circular inhomogeneity.
sigma=1/400*ones(NElement2,1);            % Make a conductivity vector.
sigma(Ind)=2/400;                         % Conductivity of the inhomogeneity.

L=16;
z=0.005*ones(L,1);
sN=max(size(Node2));
[II1,T]=Current(L,sN,'adj');

[Agrad,Kb,M,S,C]=FemMatrix(Node2,Element2,z);
A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.

[U,p,r]=ForwardSolution(NNode2,NElement2,A,C,T,[],'real'); % Simulated data.
Uel=U.Electrode(:);

Agrad1=SparseCrush(Agrad*Ind2);   % Group some of the element for the inverse computations

%%             PROCEDURE TO SOLVE THE INVERSE PROBLEM           %%

% Approximate the best homogenous resistivity.

A=UpdateFemMatrix(Agrad,Kb,M,S,ones(NElement2,1));  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);

rho0=Uref.Electrode(:)\U.Electrode(:);

A=UpdateFemMatrix(Agrad,Kb,M,S,1./rho0*ones(size(sigma)));  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);
Urefel=Uref.Electrode(:);

rho=rho0*ones(size(Agrad1,2),1);
J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField, ...
           rho,'real');

%Regularisation parameter and matrix

alpha=0.0005;
R=MakeRegmatrix(Element1);


%%%Initializations for Kalman filter and smoother
a1=10; %Coefficient for diagonal state noise covariance.
a2=0.001; %Coefficient for diagonal measurement noise covariance.
F=sparse(eye(NElement1));

%%Kalman filter and Fixed-interval smoother
[rhoF,rhoS]=Dyneit(J,Uel,L,NElement1,Urefel,a1,a2,rho0,alpha,R,F);
clf,dynadecks(rhoF,[],g1,H1,4)
figure,clf,dynadecks(rhoS,[],g1,H1,4)



