% This MATLAB program checks the feasibility of LMI from Theorem 1 of the paper
% A. Selivanov and E. Fridman, "Sampled-data relay control of semilinear diffusion PDEs," in 55th IEEE Conference on Decision and Control, 2016.

h=1e-3;         % the sampling period
q1=-5;          % lower sector bound for the nonlinearity 
q2=3;           % upper sector bound for the nonlinearity 
N=1;            % the amount of subintervals 
nu=1e-5;        % the tuning parameter (nu>0)
epsilon=1e-10;  % the shape function parameter
alpha=1;        % a desired decay rate
BC=3;           % the type of boundary conditions

K=50;           % the controller gain 
rho=0.1;        % the disturbance psize is rho*K

[M,betau,betaw]=LMI_CDC16b_th1(h,q1,q2,N,nu,epsilon,alpha,BC); 
if ~isempty(M)
    disp(['Boundary for the initial conditions: C_0=' num2str((1-rho)^2*K^2/(N*M^2))]); 
    disp(['Limit value boundary: C_\infty=' num2str((betau+betaw*rho^2)*K^2*h/(2*alpha))]); 
    disp(['The disturbance size: ' num2str(rho*K)]); 
else
    disp('LMIs are not feasible'); 
end