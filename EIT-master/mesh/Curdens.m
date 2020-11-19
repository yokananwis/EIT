function cd=Curdens(U,z,g,E);


u=U.Current;
u=u(E1);
Uel=U.Electrode;
k=1;
Uel=Uel(:,k)*ones(1,size(u,2));
cd=[zeros(size(Uel,1),2),[(Uel-u)/z(1)],zeros(size(Uel,1),2)];
cd=cd';
cd=cd(:);

%%%%%%%%%%% Mesh, do this first %%%%%%%%%%%%%%%
N=round([256, 128./r(14:-1:2),1]);
r=[14:-1:0];
eI=[8,8];
[g1,gp,H1,E1]=cirgrid_eit(r,N,eI);
Indb1=findboundary(g1,H1);

[Element2,Nodelist2]=MakeElement(H1,Indb1,E1);
[Node2]=MakeNode(Element2,Nodelist2,g1);

