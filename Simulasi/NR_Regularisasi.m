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
A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.

% This is ultimately what we want to plot:
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
lambda = 1.0e-15;
R=MakeRegmatrix(Element1);

%% First Iteration
num_iter=0;
rho=rho+(J'*J+lambda*(R'*R))\(J'*(Uel-Urefel));
rhobig=Ind2*rho;
A=UpdateFemMatrix(Agrad,Kb,M,S,1./rhobig);  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);
Urefel=Uref.Electrode(:);
J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,rho,'real');
num_iter=num_iter+1;
Obj_F=0.5*(Urefel-Uel)'*(Urefel-Uel);

figure(2)
clf,Plotinvsol(rho,g1,H1);colorbar,title(['NR Regularization Reconstruction | Iter: ' num2str(num_iter) ' steps']);drawnow;

%% Next Iterations
disp('Iterations...')
% I think the iterations just refine the guess
while Obj_F(num_iter)>0.0001
    %rho=rho+(J'*J+lambda*(R'*R))\(J'*(Uel-Urefel)-lambda*(R'*R)*rho);
    lambda=lambda/10;
    rho=rho+(J'*J+lambda*(R'*R))\(J'*(Uel-Urefel));
    rhobig=Ind2*rho;
    A=UpdateFemMatrix(Agrad,Kb,M,S,1./rhobig);  % The system matrix.
    Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);
    Urefel=Uref.Electrode(:);
    J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,rho,'real');
    num_iter=num_iter+1;
    Obj_F=[Obj_F,0.5*(Urefel-Uel)'*(Urefel-Uel)];
    
    figure(2)
    clf,Plotinvsol(rho,g1,H1);colorbar,title(['NR Regularization Reconstruction | Iter: ' num2str(num_iter) ' steps']);drawnow;
end
rho_reg=rho;

x = 1:1:numel(Obj_F);
y = Obj_F;

figure(9)
clf,plot(x,y,'-o','MarkerIndices',1:numel(y));
xticks(1:1:numel(Obj_F));
title('Fungsi Objektif NR Regularisasi');
xlabel('Iterasi ke-');
ylabel('Fungsi Objektif');