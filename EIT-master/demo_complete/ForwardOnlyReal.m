% ForwardOnlyReal demonstrates the use of 2D EIT Package for simulations with linear approximation basis
% Only the forward solution, with real impedance values.

% M. Vauhkonen 28.3.2000,
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi

load meshdata % Data for two different meshes.

NNode1=max(size(Node1));                      %The number of nodes
NElement1=max(size(Element1));                %The number of element
% The mesh for the inhomogeneities:
NNode2=max(size(Node2));                      %The number of nodes
NElement2=max(size(Element2));                %The number of elements

g1=reshape([Node1.Coordinate],2,NNode1)';
H1=reshape([Element1.Topology],3,NElement1)';
g2=reshape([Node2.Coordinate],2,NNode2)';
H2=reshape([Element2.Topology],3,NElement2)';


disp('Choose a circular inhomogeneity. Left mouse button, center, right button, radius.')
Ind1=ChooseCircle(Node2,Element2);     % Make data for an inhomogeneity.
resis = 5*ones(NElement2,1);         % Resistivity vector
resis(Ind1)=25;			           % Resistivity of the inhomogeneity.
sigma=1./resis;                        % Make a conductivity vector.
% We can use complex conductivity! I think I'd rather make it resistivity
% for ease.

figure(1)
clf,Plotinvsol(1./sigma,g2,H2);colorbar,title('Your resistivity distribution');drawnow
disp('Press any key to continue...'),pause

disp('Computes the simulated data.')
L=16;					  % The number of electrodes.
z=0.005*ones(L,1);			  % Contact impedances.
% Specify the current pattern. Set the last argument to L/2 for 'opp', L-1
% for 'tri', or no argument. 
% The current pattern should be 'tri' to get correct magnitude results.
% We want the rms current to be 800uA.
rms = 800e-6;

[II1,T]=Current(L,NNode2,'tri',rms);	  

[Agrad,Kb,M,S,C]=FemMatrix(Node2,Element2,z);
A=UpdateFemMatrix(Agrad,Kb,M,S,sigma);  % The system matrix.

[U,p,r]=ForwardSolution(NNode2,NElement2,A,C,T,[],'real'); % Simulated data.
Uel=U.Electrode(:);
figure(2)
clf,Plotinvsol(U.Current,g2,H2); colorbar; drawnow;

