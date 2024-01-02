%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%   Exercise 14: Adjustment Calculation - part IX  
% 
%   Author         : Anastasia Pasioti
%   Version        : January 29, 2018
%   Last changes   : February 14, 2022
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;
format long g;

%--------------------------------------------------------------------------
%   Task 1
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Vector of observations
L = [ -4.0; 8.0; 7.7];

%Number of observations
no_n = length(L);

%Initial values for the unknowns
x = sqrt(15.7/4);
y = sqrt(8.0 - 15.7/4);

%Vector of initial values for the unknowns
X_0 = [x y]';

%Number of unknowns
no_u = length(X_0);

%Number of constraints
no_b = 1;

%Redundancy
r = no_n - no_u + no_b;   

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
%S_LL =    

%Theoretical standard deviation
% sigma_0 =      %a priori

%Cofactor matrix of the observations
Q_LL = eye(no_n);

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-7;
delta = 10^-12;
max_x_hat = 10^Inf;

%Number of iterations
iteration = 0;

while max_x_hat>epsilon || Check2>delta            
    
     %Observations as functions of the approximations for the unknowns
     L_0(1) =x + y - 2*y^2;
     L_0(2) = x^2 + y^2;
     L_0(3) = 3*x^2-y^2;

     %Vector of reduced observations
     l = L-L_0';
    
     %Design matrix A with the elements from the Jacobian matrix J
     A(1,1) = 1;
     A(1,2) = 1-4*X_0(2);

     A(2,1) = 2*X_0(1);
     A(2,2) = 2*X_0(2);

     A(3,1) = 6*X_0(1);
     A(3,2) = -2*X_0(2);
     
     %Design matrix C with the elements from the Jacobian matrix J
    %C = [3*X_0(1)^2 (X_0(2)^(-2/3))/3]; %1st constraint
     C = [-1/(X_0(1)^2) 2*(X_0(2))];     %2st constraint

     %Normal matrix
     N = A'*P*A;
     
     %Extended normal matrix
     N_ext = [N C'
              C 0];
     
     %Vector of right hand side of normal equations
     n = A'*P*l;
     
     %Extended vector of right hand side of normal equations
     %n_ext = [n; 9.0-X_0(1)^3-X_0(2)^(1/3)]; %1st constraint
     n_ext = [         n
              5.0-(1/X_0(1))-X_0(2)^2]; %2st constraint

     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx_ext = inv(N_ext);
     Q_xx = Q_xx_ext(1:no_u,1:no_u);

     %Solution of the normal equations
     x_hat = Q_xx_ext*n_ext;
       
     %Update
     X_0 = X_0 + x_hat(1:no_u);
    
     x = X_0(1);
     y = X_0(2);
     
     %Lagrange multiplier
     k = x_hat(end);
    
     %Check 1
     max_x_hat = max(abs(x_hat(1:no_u)));
     
     %Vector of residuals
     v = A*x_hat(1:no_u) - l;
 
     %Vector of adjusted observations
     L_hat = L+v;
    
     %Function
     vTPv = v'*P*v;
    
     %Functional relationships without the observations
     phi_X_hat = [x+y-2*y^2; 
                  x^2+y^2; 
                  3*x^2-y^2];

     %Check 2
     Check2 = max(abs(L_hat-phi_X_hat));
    
     %Update number of iterations
     iteration = iteration+1;
  
end

if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end


x
y
%x^3+y^(1/3)
(1/x)+y^2

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));        

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));     

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));         

    

