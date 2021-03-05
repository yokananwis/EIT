% EidorsDemo1 Demonstrates the use of 2D EIT Package for simulations with linear approximation basis

% M. Vauhkonen 28.3.2000,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

clear;
load meshdata % Data for two different meshes.

NNode1=max(size(Node1));                      %The number of nodes
NElement1=max(size(Element1));                %The number of element
% I think this is for the single inhomogeneity:
NNode2=max(size(Node2));                      %The number of nodes
NElement2=max(size(Element2));                %The number of elements

g1=reshape([Node1.Coordinate],2,NNode1)';
H1=reshape([Element1.Topology],3,NElement1)';
g2=reshape([Node2.Coordinate],2,NNode2)';
H2=reshape([Element2.Topology],3,NElement2)';


disp('Generating Finite Element Model.')
Ind=ChooseCircle(Node2,Element2);       % Make data for an inhomogeneity.
sigma=1/400*ones(NElement2,1);            % Make a conductivity vector.
sigma(Ind)=1/200;			  % Conductivity of the inhomogeneity.
% sigma = CreateInhomogeneities(Node2,Element2,7);

% Eventually we'll want to get rid of Plotinvsol or rewrite it.
figure(1)
clf,Plotinvsol(1./sigma,g2,H2);colorbar,title('Finite Element Model');
disp('Press any key to continue...'),pause

disp('Computes the simulated data.')
L=16;					  % The number of electrodes.
z=0.005*ones(L,1);			  % Contact impedances.
rms = 800e-6;
[II1,T]=Current(L,NNode2,'tri',rms);	  % Trigonometric current pattern.

[Agrad,Kb,M,S,C]=FemMatrix(Node2,Element2,z);
Sys_Mat=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.

% This is ultimately what we want to plot:
[U,p,r]=ForwardSolution(NNode2,NElement2,Sys_Mat,C,T,[],'real'); % Simulated data.
Uel=U.Electrode(:);

Agrad1=Agrad*Ind2;   % Group some of the element for the inverse computations


%%             PROCEDURE TO SOLVE THE INVERSE PROBLEM           %%

% Approximate the best homogenous resistivity.


disp('Solves the full nonlinear inverse problem by Landweber iteration.')

disp('Initialisations...')

Sys_Mat=UpdateFemMatrix(Agrad,Kb,M,S,ones(NElement2,1));  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,Sys_Mat,C,T,[],'real',p,r);

rho0=Uref.Electrode(:)\U.Electrode(:);

Sys_Mat=UpdateFemMatrix(Agrad,Kb,M,S,1./rho0*ones(size(sigma)));  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,Sys_Mat,C,T,[],'real',p,r);
Urefel=Uref.Electrode(:);

% initial resistivity distribution
rho=rho0*ones(size(Agrad1,2),1);

% initial estimation
J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,rho,'real');
max_eigen=eigs(J'*J,1);
alpha=2/max_eigen;
p_norm=norm(alpha*(J'*J),2);
num_iter=0;

disp('Iterations...')
% I think the iterations just refine the guess
while p_norm >= 2
 %Calculate Jacobian (Sensitivity Matrix)
 J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,rho,'real');
 rhobig=Ind2*rho;
 Sys_Mat=UpdateFemMatrix(Agrad,Kb,M,S,1./rhobig);  % The system matrix.
 Uref=ForwardSolution(NNode2,NElement2,Sys_Mat,C,T,[],'real',p,r);
 Urefel=Uref.Electrode(:);
 max_eigen=eigs(J'*J,1);
 alpha=2/max_eigen;
 p_norm=norm(alpha*(J'*J),2);
 rho=rho+alpha*J'*(Uel-Urefel);
 %rho=rho+alpha*J'*(Uel-J*rho);
 
 % resistivity distribution plot
 figure(5)
 clf,Plotinvsol(rho,g1,H1);colorbar,title(['Iterasi Landweber | Iter: ' num2str(num_iter) ' steps']);drawnow;
 num_iter = num_iter + 1;
end