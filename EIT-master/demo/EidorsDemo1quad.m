%EidorsDemo1quad Demonstrates the use of 2D EIT Package with quadratic basis
% EidorsDemo1quad Demonstrates the use of 2D EIT Package for simulations with quadratic basis.
%
% M. Vauhkonen 28.3.2000,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

load meshdata2nd % Data for two different meshes.

NNode1=max(size(Node));                      %The number of nodes
NElement1=max(size(Element));                %The number of element
NNode2=max(size(Node1));                      %The number of nodes
NElement2=max(size(Element1));                %The number of elements

g1=reshape([Node.Coordinate],2,NNode1)';
H1=reshape([Element.Topology],3,NElement1)';
g2=reshape([Node1.Coordinate],2,NNode2)';
H2=reshape([Element1.Topology],6,NElement2)';


Ind=ChooseElements(Node,Element,[],10);   % Make data for an inhomogeneity.
sigma=1/400*ones(NElement2,1);            % Make a conductivity vector.
sigma(Ind)=2/400;			  % Conductivity of the inhomogeneity.
L=16;					  % The number of electrodes.
z=0.005*ones(L,1);			  % Contact impedances.
[II1,T]=Current(L,NNode2,'tri');	  % Trigonometric current pattern.

[Agrad,Kb,M,S,C]=FemMatrix(Node1,Element1,z);
A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.

[U,p,r]=ForwardSolution(NNode2,NElement2,A,C,T,[],'real'); % Simulated data.
Uel=U.Electrode(:);


%%             PROCEDURE TO SOLVE THE INVERSE PROBLEM           %%

% Approximate the best homogenous resistivity.

A=UpdateFemMatrix(Agrad,Kb,M,S,ones(NElement2,1));  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);

rho0=Uref.Electrode(:)\U.Electrode(:);
rho=rho0*ones(size(sigma));

A=UpdateFemMatrix(Agrad,Kb,M,S,1./rho);  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);
Urefel=Uref.Electrode(:);

J=Jacobian(Node1,Element1,Agrad,Uref.Current,Uref.MeasField,rho,'real');

% Regularisation parameter and matrix

alpha=0.000005; 
R=MakeRegmatrix(Element);

iter=5;

for ii=1:iter
 rho=rho+(J'*J+alpha*R'*R)\(J'*(Uel-Urefel)-alpha*R'*R*rho);
 A=UpdateFemMatrix(Agrad,Kb,M,S,1./rho);  % The system matrix.
 Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);
 Urefel=Uref.Electrode(:);
 J=Jacobian(Node1,Element1,Agrad,Uref.Current,Uref.MeasField,rho,'real');
 clf,Plotinvsol(rho,g1,H1);drawnow;
end






