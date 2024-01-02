%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%    Exercise 7: Adjustment Calculation - part II  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 09, 2018
%   Last changes   : December 14, 2022
%
%--------------------------------------------------------------------------

clc;
close all;
format long;

%--------------------------------------------------------------------------
%   Task 1: Adjustment of a straight line
%--------------------------------------------------------------------------
%  equation of a line 

%  y = mx + c --------------- functional model

%--------------------------------------------------------------------------
%   Observations and redundancy
%--------------------------------------------------------------------------
%Vector of observations 
L = [0.1 1.1 1.8 2.4 ]';
%Number of observations
no_n = 4; 
%Number of unknowns
no_u = 2;
%Redundancy
r = no_n - no_u;
%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%V-C Matrix of the observations
S_ll = ([0.02^2;0.01^2;0.04^2;0.02^2].* eye(4));
%Theoretical reference standard deviation
sigma_0 = 0.001;      %a priori

%Cofactor matrix of the observations
Q_ll = ((sigma_0^2)^-1 * S_ll);
%Weight matrix
P = Q_ll^-1;
%--------------------------------------------------------------------------
%  Adjustment ( A'PAX =A'PL )
%--------------------------------------------------------------------------
%Design matrix
A =[1 1;
    2 1;
    3 1;
    4 1]; 
%Normal matrix
N = A'*P*A;       
%Vector of right hand side of normal equations
n = A'*P*L;
%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_xx = N^-1;
%Solution of normal equation
X_hat = Q_xx*n;
%Estimated unknown parameters
m =X_hat(1);
c =X_hat(2);
%Vector of residuals
v = (A*X_hat)-L ;
%Objective function
vTPv = v'*P*v; 
%Vector of adjusted observations
L_hat = L+v;
%Final check
if L_hat == (A*X_hat)
    disp('No errrors found!')
else
    disp('crosscheck your steps again!')
end

%Empirical reference standard deviation
s_0 = sqrt(vTPv/(no_n-no_u));       %a posteriori
%VC matrix of adjusted unknowns
E_xx_hat =(s_0^2*Q_xx);
%Standard deviation of the adjusted unknowns
s_m = sqrt(E_xx_hat(1));
s_c = sqrt(E_xx_hat(2,2));
%Cofactor matrix of adjusted observations
Q_ll_hat = A*Q_xx*A';
%VC matrix of adjusted observations
E_ll_hat1 = (s_0^2*Q_ll_hat);
%Standard deviation of the adjusted observations
s_L_hat1 = sqrt(E_ll_hat1(1));  %-----------obs1
s_L_hat2 = sqrt(E_ll_hat1(2,2));%-----------obs2
s_L_hat3 = sqrt(E_ll_hat1(3,3));%-----------obs3
s_L_hat4 = sqrt(E_ll_hat1(4,4));%-----------obs4
%Cofactor matrix of the residuals
Q_vv = Q_ll - Q_ll_hat;
%VC matrix of residuals
E_vv = (s_0^2*Q_vv);
%Standard deviation of the residuals
s_v1 = sqrt(E_vv(1));
s_v2 = sqrt(E_vv(2,2));
s_v3 = sqrt(E_vv(3,3));
s_v4 = sqrt(E_vv(4,4));












