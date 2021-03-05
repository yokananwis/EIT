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
Uel=xlsread('E16_Center.xlsx');
Uel=Uel(:);

%%             PROCEDURE TO SOLVE THE INVERSE PROBLEM           %%

% Approximate the best homogenous resistivity.

MeasPatt=-T;
[Agrad,Kb,M,S,C]=FemMatrix(Node2,Element2,z);
Agrad1=Agrad*Ind2;
A=UpdateFemMatrix(Agrad,Kb,M,S,ones(NElement2,1));  % The system matrix.
[Uref,p,r]=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real');
Urefel=Uref.Electrode(:);
Urefel=RemoveVolt(Urefel,L);

rho0=Urefel(:)\Uel(:);

A=UpdateFemMatrix(Agrad,Kb,M,S,1/rho0*ones(NElement2,1));  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real',p,r);
Urefel=Uref.Electrode(:);
Urefel=RemoveVolt(Urefel,L);
Urefel=Urefel(:);

rho=rho0*ones(size(Ind2,2),1);
J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField, ...
           rho,'real');
J=RemoveJacob(J,L);

%Regularisation parameter and matrix
lambda = 0.001;
R=MakeRegmatrix(Element1);

%% First Iteration
num_iter=0;

rho=rho+(J'*J+lambda*(R'*R))\(J'*(Uel-Urefel));
rhobig=Ind2*rho;
A=UpdateFemMatrix(Agrad,Kb,M,S,1./rhobig);  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real',p,r);
Urefel=Uref.Electrode(:);
Urefel=RemoveVolt(Urefel,L);
Urefel=Urefel(:);
J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,rho,'real');
J=RemoveJacob(J,L);
num_iter=num_iter+1;
Obj_F=0.5*(Urefel-Uel)'*(Urefel-Uel);

figure(1)
clf,Plotinvsol(rho,g1,H1);colorbar,title(['NR Regularization Reconstruction | Iter: ' num2str(num_iter) ' steps']);drawnow;

%% Next Iterations
disp('Iterations...')
% I think the iterations just refine the guess
while Obj_F(num_iter)>0.01
    lambda=lambda/10;
    rho=rho+(J'*J+lambda*(R'*R))\(J'*(Uel-Urefel));
    rhobig=Ind2*rho;
    A=UpdateFemMatrix(Agrad,Kb,M,S,1./rhobig);  % The system matrix.
    Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'real',p,r);
    Urefel=Uref.Electrode(:);
    Urefel=RemoveVolt(Urefel,L);
    Urefel=Urefel(:);
    J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,rho,'real');
    J=RemoveJacob(J,L);
    num_iter=num_iter+1;
    Obj_F=[Obj_F,0.5*(Urefel-Uel)'*(Urefel-Uel)];
    figure(1)
    clf,Plotinvsol(rho,g1,H1);colorbar,title(['NR Regularization Reconstruction | Iter: ' num2str(num_iter) ' steps']);drawnow;
end
rho_reg=rho;

x = 1:1:numel(Obj_F);
y = Obj_F;

figure(2)
clf,plot(x,y,'-o','MarkerIndices',1:numel(y));
xticks(1:1:numel(Obj_F));
title('Fungsi Objektif NR Regularisasi');
xlabel('Iterasi ke-');
ylabel('Fungsi Objektif');