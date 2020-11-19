%EidorsDemoCompl Demonstrates the use of 2D EIT Package for reconstructing admittivity
% EidorsDemoCompl Demonstrates the use of 2D EIT Package for reconstructing admittivity

% M. Vauhkonen 20.5.2000,
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


Ind=ChooseCircle(Node2,Element2); % Make data for an inhomogeneity.
sigma=1/400*ones(NElement2,1);            % Make a conductivity vector.
sigma(Ind)=2/400;			  % Conductivity of the inhomogeneity.

Ind=ChooseCircle(Node2,Element2); % Make data for an inhomogeneity.
we=1/400*ones(NElement2,1);               % Make a scaled permittivity vector.
we(Ind)=0.5/400;                            % Conductivity of the inhomogeneity.

adm=sigma+i*we;                           % Admittivity. 

L=16;					  % The number of electrodes.
z=0.005*ones(L,1);			  % Contact impedances.
[II1,T]=Current(L,NNode2,'adj');	  % Adjacent current pattern.

[Agrad,Kb,M,S,C]=FemMatrix(Node2,Element2,z);
A=UpdateFemMatrix(Agrad,Kb,M,S,adm);  % The system matrix.

U=ForwardSolution(NNode2,NElement2,A,C,T,[],'comp'); % Simulated data.
Uel=U.Electrode(:);

Agrad1=SparseCrush(Agrad*Ind2);   % Group some of the elements for the inverse computations


%%             PROCEDURE TO SOLVE THE INVERSE PROBLEM           %%

% Start from the background value

adm=1/400*ones(NElement1,1)+i*1/400*ones(NElement1,1);
admbig=Ind2*adm;
A=UpdateFemMatrix(Agrad,Kb,M,S,admbig);  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'comp');
Urefel=Uref.Electrode(:);
J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,[],'comp');
adm=[real(adm);imag(adm)];

%Regularisation parameter and matrix

alpha=0.9; 
R=MakeRegmatrix(Element1);
R=[R,zeros(size(R));zeros(size(R)),R];

iter=5;

for ii=1:iter
 adm=adm+(J'*J+alpha*R'*R)\(J'*(Uel-Urefel)-alpha*R'*R*adm);
 adm=reshape(adm,max(size(Element1)),2);
 adm=adm(:,1)+i*adm(:,2);
 admbig=Ind2*adm;
 A=UpdateFemMatrix(Agrad,Kb,M,S,admbig);  % The system matrix.
 Uref=ForwardSolution(NNode2,NElement2,A,C,T,[],'comp');
 Urefel=Uref.Electrode(:);
 J=Jacobian(Node2,Element2,Agrad1,Uref.Current,Uref.MeasField,[],'comp');
 figure(1),clf,Plotinvsol(real(adm),g1,H1),title('Conductivity');colorbar,drawnow
 figure(2),clf,Plotinvsol(imag(adm),g1,H1),title('Scaled permittivity');colorbar,drawnow,
 adm=[real(adm);imag(adm)];ii
end






