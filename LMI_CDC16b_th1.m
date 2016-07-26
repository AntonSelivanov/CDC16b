function [M,betau,betaw]=LMI_CDC16b_th1(h,q1,q2,N,nu,epsilon,alpha,BC)
% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper 
% A. Selivanov and E. Fridman, "Sampled-data relay control of semilinear diffusion PDEs," in 55th IEEE Conference on Decision and Control, 2016.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% h       - the sampling period; 
% q1,q2   - the sector bound for the nonlinearity; 
% N       - the amount of subintervals; 
% nu      - the tuning parameter (nu>0); 
% epsilon - the shape function parameter; 
% alpha   - a desired decay rate; 
% BC      - the type of boundary conditions: 1 - Dirichlet, 2 - Neumann, 3 - both. 

% Output: 
% M, betau, betaw - the parameters required to calculate the initial and limit bounds. 

%% Decision variables 
p=sdpvar; 
q=sdpvar; 
Mvar=sdpvar; 
lambdaphi=sdpvar;
lambdakappa=sdpvar; 
if BC==2
    lambda=0; 
else
    lambda=sdpvar; 
end
betauvar=sdpvar; 
betawvar=sdpvar; 

%% Notations
if BC==1
    d=4; 
else
    d=1;
end

%% The main LMI (Xi<=0)
Xi=blkvar; 
Xi(1,1)=-2*Mvar+2*lambdakappa*epsilon*N^3*pi^2/nu+2*alpha-lambdaphi*q1*q2-d*lambda*pi^2/4; 
Xi(1,2)=1+lambdaphi*(q1+q2)/2; 
Xi(1,3)=Mvar; 
Xi(1,5)=h*Mvar; 
Xi(2,2)=-lambdaphi; 
Xi(2,4)=-q*h;
Xi(2,8)=p*h*exp(alpha*h); 
Xi(3,3)=-lambdakappa*(N*pi)^2/(1+nu); 
Xi(3,5)=-h*Mvar; 
Xi(4,4)=-2*q*h; 
Xi(4,6)=-q*h; 
Xi(4,7)=-q*h; 
Xi(4,8)=p*h*exp(alpha*h); 
Xi(5,5)=-p*h*pi^2/4; 
Xi(5,6)=h; 
Xi(5,7)=h; 
Xi(6,6)=-betauvar*h; 
Xi(6,8)=p*h*exp(alpha*h); 
Xi(7,7)=-betawvar*h; 
Xi(7,8)=p*h*exp(alpha*h); 
Xi(8,8)=-p*h; 
Xi=sdpvar(Xi); 

%% Solution of LMIs
LMIs=[Xi<=0, 2*alpha*q*h+lambdakappa+lambda<=2, p>=0, q>=0, Mvar>=0, lambdaphi>=0, lambdakappa>=0, lambda>=0, betauvar>=0, betawvar>=0]; 
options=sdpsettings('solver','lmilab','verbose',0);
sol=optimize(LMIs,betauvar+betawvar*.1^2,options); 

M=[]; betau=[]; betaw=[]; 
if sol.problem == 0
    [primal,~]=check(LMIs); % Checking that the solver returned a proper solution
    if min(primal)>=0 
        M=double(Mvar); 
        betau=double(betauvar); 
        betaw=double(betawvar); 
    end
else
    yalmiperror(sol.problem); 
end