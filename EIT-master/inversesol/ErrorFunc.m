function no=ErrorFunc(rho,Uel,Agrad,Kb,M,S,R,alpha,NNode2,NElement2,C,T,MeasPatt,p,r,L,Ind2);

A=UpdateFemMatrix(Agrad,Kb,M,S,Ind2*rho);  % The system matrix.
Uref=ForwardSolution(NNode2,NElement2,A,C,T,MeasPatt,'real',p,r);
Urefel=Uref.Electrode(:);
Urefel=RemoveVolt(Urefel,L);
no=norm((Uel(:)-Urefel(:)))^2+alpha*norm(R*rho)^2;

