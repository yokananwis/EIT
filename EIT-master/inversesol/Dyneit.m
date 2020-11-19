function [rhoF,rhoS]=Dyneit(J,Uelmit,L,nump,Uellask,a1,a2,rho0,alpha,R,F);

%Dyneit Estimates resistivity distributions by Kalman filter and Kalman smoother 
% Function [rhoF,rhoS]=Dyneit(J,Uelmit,L,nump,Uellask,a1,a2,rho0,alpha,R,F);
% estimates resisitivity distributions by Kalman filter and smoother.
%
% INPUT
%
% J = Jacobian
% Uelmit = measured voltages on the electrodes
% L = the number of electrodes
% nump = the number of parameters
% Uellask = computed voltages on the electrodes at the linearisation point 
% a1 = state noise covariance coefficient
% a2 = measurement noise covariance coefficient
% rho0 = the first guess for the resistivity
% alpha = regularization parameter
% R = regularization matrix
% F = state transition matrix
%
% OUTPUT
%
% rhoF = resistivity estimate computed by the Kalman filter
% rhoS = resistivity estimate computed by the Kalman smoother 

% P.J. Vauhkonen 5.8.1999
% Modified for EIDORS by M. Vauhkonen 17.11.1999
% University of Kuopio, Department of Applied Physics, PO Box 1627,
% FIN-70211 Kuopio, Finland, email: Marko.Vauhkonen@uku.fi


Cw=a1*sparse(eye(nump));
Cv1=[a2*sparse(eye(L)),sparse(zeros(L,nump))];
Cv2=[sparse(zeros(nump,L)),sparse(eye(nump))];
Cv=[Cv1;Cv2]; 
rhoP=zeros(nump,L);
rhoF=zeros(nump,L);
A=zeros(nump,nump*(L-1));
rhoP(:,1)=rho0;
CP=sparse(eye(nump));

%%%%%% Filter

for T=1:L
 z1=[Uelmit((T-1)*L+1:T*L)-Uellask((T-1)*L+1:T*L);zeros(nump,1)];
 J2=J((T-1)*L+1:T*L,:);
 J1=[J2;alpha*R];
 K=CP*J1'*inv(J1*CP*J1'+Cv);
 CF=(eye(nump,nump)-K*J1)*CP;
 rhoF(:,T)=rhoP(:,T)+K*(z1-J1*(rhoP(:,T)-rho0));
 if T<L
  CP=F*CF*F'+Cw;
  rhoP(:,T+1)=F*rhoF(:,T);
  A(:,(T-1)*nump+1:T*nump)=CF*F'*inv(CP);
 else
 end
T
end

%%%%% Smoother
rhoS=zeros(nump,L);
rhoS(:,L)=rhoF(:,L);

for t=L:-1:2
rhoS(:,t-1)=rhoF(:,t-1)+A(:,(t-2)*nump+1:(t-1)*nump)*(rhoS(:,t)-rhoP(:,t));
end


















