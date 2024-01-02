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
%   Task 2
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and redundancy
%--------------------------------------------------------------------------
%Height differences/observations (<h+v =Hb-Ha)
a1=10.509;
a2=5.360;
a3=-8.523;
a4=-7.348;
a5=-3.167;
a6=15.881;
L =[a1 a2 a3 a4 a5 a6]';      %[m]

%Benchmarks - error free
Ha = 123.456;         %[m]
Hb = 654.321;         %[m]

%Vector of observations with benchmarks
I_a =[(a1+Ha)  a2    a3 (a4-Ha)  a5    (a6+Ha)]'; % with benchmark A
I_b =[(a1-Hb) a2+Hb  a3   a4   (a5+Hb)    a6]';   % with benchmark B

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
S_LL = eye(6).*[0.006^2  0.004^2  0.005^2  0.003^2  0.004^2  0.012^2];

%Theoretical standard deviation
sigma_0 = 1;  %a priori

%Cofactor matrix of the observations
Q_LL = sigma_0^-2 * S_LL;  
%Weight matrix
P = Q_LL^-1; 

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Design matrix
A_Ha =[1 0 0;-1  1  0;0 -1  1;0  0 -1;-1  0  1;0  1  0];
A_Hb =[-1 0 0; 0  1  0;0 -1  1;1  0 -1; 0  0  1;-1 1  0];
    
%Normal matrix
N_Ha =A_Ha'*P*A_Ha; % for benchmark A
N_Hb =A_Hb'*P*A_Hb; % for benchmark B       
%Vector of right hand side of normal equations
n_Ha =A_Ha'*P*I_a; % for beachmark A
n_Hb =A_Hb'*P*I_b; % for benchmark B

%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_XX_Ha =N_Ha^-1;  % for benchmark A
Q_XX_Hb =N_Hb^-1;  % for benchmark B
%Solution of normal equation
X_hat_Ha =Q_XX_Ha*n_Ha; % for benchmark A
X_hat_Hb =Q_XX_Hb*n_Hb; % for benchmark B  
%Estimated unknown parameters
Hb_a =X_hat_Ha(1);     %for point A with beachmark A 
Hc_a =X_hat_Ha(2);     %for point B with benchmark A
Hd_a =X_hat_Ha(3);     %for point C with benchmark A
Ha_b =X_hat_Hb(1);     %for point A with benchmark B
Hc_b =X_hat_Hb(2);     %for point C with benchmark B
Hd_b =X_hat_Hb(3);     %for point D with benchmark B
%Vector of residuals
v_Ha = A_Ha*X_hat_Ha-I_a; % with benchmark A
v_Hb = A_Hb*X_hat_Hb-I_b; % with benchmark B
%Objective function
vTPv_Ha =v_Ha'*P*v_Ha; % with benchmark A
vTPv_Hb =v_Hb'*P*v_Hb; % with benchmark B

%Vector of adjusted observations

L_hat1 =L+v_Ha; % with benchmark A
L_hat2 =L+v_Hb; % with benchmark B

%Final check (to check for computational errors)

L1_hat_Ha =I_a+v_Ha;
L1_hat_Hb =I_b+v_Hb;
Phi_X_hat_Ha =A_Ha*X_hat_Ha; 
Phi_X_hat_Hb =A_Hb*X_hat_Hb;

if all(L1_hat_Ha - Phi_X_hat_Ha < 10^-5) ||all( L1_hat_Hb - Phi_X_hat_Hb < 10^-5)
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

%Empirical reference standard deviation
s_0_Ha = sqrt(vTPv_Ha/r);      %a posteriori
s_0_Hb = sqrt(vTPv_Hb/r);      %a posteriori

%VC matrix of adjusted unknowns
S_XX_hat_Ha = s_0_Ha^2*Q_XX_Ha;
S_XX_hat_Hb = s_0_Hb^2*Q_XX_Hb;
%Standard deviation of the adjusted unknowns

s_X_Ha =sqrt(diag(S_XX_hat_Ha));
s_X_Hb =sqrt(diag(S_XX_hat_Hb));

%Cofactor matrix of adjusted observations

Q_LL_hat_Ha = A_Ha*Q_XX_Ha*A_Ha';  % for benchmark A
Q_LL_hat_Hb = A_Hb*Q_XX_Hb*A_Hb';  % for benchmark B

%VC matrix of adjusted observations
S_LL_hat_Ha = s_0_Ha^2*Q_LL_hat_Ha; 
S_LL_hat_Hb = s_0_Hb^2*Q_LL_hat_Hb; 

%Standard deviation of the adjusted observations

s_L_hat_Ha =sqrt(diag(S_LL_hat_Ha));  % for benchmark A
s_L_hat_Hb =sqrt(diag(S_LL_hat_Hb));  % for benchmark B

%Cofactor matrix of the residuals

Q_vv_Ha = Q_LL - Q_LL_hat_Ha;   % for benchmark A
Q_vv_Hb = Q_LL - Q_LL_hat_Hb;   % for benchmark B

%VC matrix of residuals

S_vv_Ha = s_0_Ha^2 * Q_vv_Ha;   % for benchmark A
S_vv_Hb = s_0_Hb^2 * Q_vv_Hb;   % for benchmark B

%Standard deviation of the residuals
s_v_Ha =sqrt(diag(S_vv_Ha)); 
s_v_Hb =sqrt(diag(S_vv_Hb)); 

