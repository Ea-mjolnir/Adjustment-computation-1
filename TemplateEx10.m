%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%  Exercise 10: Adjustment Calculation - part V  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 11, 2018
%   Last changes   : January 18, 2023
%
%--------------------------------------------------------------------------

clc;
close all;

%--------------------------------------------------------------------------
%   Task 1
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and redundancy
%--------------------------------------------------------------------------
%Height differences/observations (<H->ab =Hb-Ha)
a1=5.1;
a2=2.340;
a3=-1.250;
a4=-6.130;
a5=-0.68;
a6=-3.0;
a7=1.70;
L =[a1 a2 a3 a4 a5 a6 a7]';      %[m]

%Benchmarks - error free
H100 = 100;            %[m]
H200 = 107.50;         %[m]

%Vector of observations with benchmarks
I =[(a1+100) (a2-107.5) (a3+107.5) (a4-100) a5 (a6+107.5) a7]'; 

%Number of observations
no_n = length(L);

%Number of unknowns
no_u =3;

%Redundancy
r = no_n - no_u; 

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
S_LL = eye(7);

%Theoretical standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_LL = sigma_0^-2 * S_LL;

%Weight matrix
P = Q_LL^-1; 

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Design matrix
A =[1  0  0
   -1  0  0
    0  0  1
    0  0 -1
   -1  1  0
    0  1  0
    0 -1  1];
    
%Normal matrix
N =A'*P*A; 
        
%Vector of right hand side of normal equations
n =A'*P*I;
    
%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_XX =N^-1;
    
%Solution of normal equation
X_hat =Q_XX*n;
    
%Estimated unknown parameters
HA =X_hat(1);     %for point A  
HB =X_hat(2);     %for point B
HC =X_hat(3);     %for point C
     
%Vector of residuals
v = A*X_hat-I;

%Objective function
vTPv =v'*P*v; 

%Vector of adjusted observations
L_hat =L+v; 

%Final check -- (to check for computational errors)
L1_hat =I+v;
Phi_X_hat = A*X_hat; 

if L1_hat-Phi_X_hat<10^-5
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);      %a posteriori

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_XX;

%Standard deviation of the adjusted unknowns
s_X =sqrt(diag(S_XX_hat));

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_XX*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat; 

%Standard deviation of the adjusted observations
s_L_hat =sqrt(diag(S_LL_hat));

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v =sqrt(diag(S_vv)); 

disp('Observation : L  L_hat  s_L_hat  v  s_v ')
disp('Unknowns : HA HB HC')
unknowns=[HA s_X(1)
          HB s_X(2)
          HC s_X(3)]
obs = [L L_hat s_L_hat v s_v]




