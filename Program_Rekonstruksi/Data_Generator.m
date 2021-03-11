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

L=16;					  % The number of electrodes.
z=0.005*ones(L,1);			  % Contact impedances.
[II1,T]=Current(L,NNode2,'adj');	  % Adjacent current pattern.

disp('Generating Finite Element Model.')
Ind=DrawCircle(Node2,Element2);       % Make data for an inhomogeneity.
sigma=1/400*ones(NElement2,1);            % Make a conductivity vector.
sigma(Ind)=1/200;			  % Conductivity of the inhomogeneity.
% sigma = CreateInhomogeneities(Node2,Element2,7);

% Eventually we'll want to get rid of Plotinvsol or rewrite it.
figure(3); clf;
clf,Plotinvsol(1./sigma,g2,H2);colorbar,title('Your resistivity distribution');

disp('Computes the simulated data.')
L=16;					  % The number of electrodes.
z=0.005*ones(L,1);			  % Contact impedances.
rms = 800e-6;
[II1,T]=Current(L,NNode2,'tri',rms);	  % Trigonometric current pattern.

load adj_meas_pattern
MeasPatt=-T;
[Agrad,Kb,M,S,C]=FemMatrix(Node2,Element2,z);
A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.

% This is ultimately what we want to plot:
[U,p,r]=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real'); % Simulated data.
Uel=U.Electrode;
Uel=RemoveVolt(Uel,L);

filename='E16_Generated_Data.xlsx';
delete(filename);
xlswrite(filename, Uel);
