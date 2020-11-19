%EidorsDemo1 Demonstrates the use of 2D EIT Package with linear basis
% EidorsDemo1 Demonstrates the use of 2D EIT Package for simulations with linear approksimation basis

% M. Vauhkonen 28.3.2000,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

load "D:/Kuliah S2/Tesis Tercinta/Simulasi/EIT-master/data/meshdata.mat" % Data for two different meshes.

NNode1=max(size(Node1));                      %The number of nodes
NElement1=max(size(Element1));                %The number of element
NNode2=max(size(Node2));                      %The number of nodes
NElement2=max(size(Element2));                %The number of elements

g1=reshape([Node1.Coordinate],2,NNode1)';
H1=reshape([Element1.Topology],3,NElement1)';
g2=reshape([Node2.Coordinate],2,NNode2)';
H2=reshape([Element2.Topology],3,NElement2)';


disp('Choose a circular inhomogeneity. Left mouse button, center, right button, radius.')
Ind=ChooseCircle(Node2,Element2);       % Make data for an inhomogeneity.
sigma=1/400*ones(NElement2,1);            % Make a conductivity vector.
sigma(Ind)=2/400;			  % Conductivity of the inhomogeneity.

clf,Plotinvsol(1./sigma,g2,H2);colorbar,title('Your resistivity distribution');drawnow
disp('Press any key to continue...'),pause

disp('Computes the simulated data.')
L=16;					  % The number of electrodes.
z=0.005*ones(L,1);			  % Contact impedances.
[II1,T]=Current(L,NNode2,'tri');	  % Trigonometric current pattern.

[Agrad,Kb,M,S,C]=FemMatrix(Node2,Element2,z);
A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.

[U,p,r]=ForwardSolution(NNode2,NElement2,A,C,T,[],'real'); % Simulated data.
Uel=U.Electrode(:);

Agrad1=Agrad*Ind2;   % Group some of the element for the inverse computations


%%             PROCEDURE TO SOLVE THE INVERSE PROBLEM           %%

% Approximate the best homogenous resistivity.


disp('Solves the full nonlinear inverse problem by regularised Gauss-Newton iteration.')

disp('Initialisations...')

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

alpha = 0.000005; 
R=MakeRegmatrix(Element1);

iter=5;

disp('Iterations...')

for ii=1:iter
 rho=rho+(J'*J+alpha*R'*R)\(J'*(Uel-Urefel)-alpha*R'*R*rho);
 rhobig=Ind2*rho;
 A=UpdateFemMatrix(Agrad,Kb,M,S,1./rhobig);  % The system matrix.
 Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);
 Urefel=Uref.Electrode(:);
 J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,rho,'real');
 clf,Plotinvsol(rho,g1,H1);colorbar,title([num2str(ii) '. step']);drawnow;
end










