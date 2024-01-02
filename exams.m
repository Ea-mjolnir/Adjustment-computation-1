%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%    Exercise 9: Adjustment Calculation - part IV  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 11, 2018
%   Last changes   : January 11, 2023
%
%--------------------------------------------------------------------------

clc;
close all;

%--------------------------------------------------------------------------
%   a+(b*x)+(c*x.^2)-sin(x)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Load data
data = load('data.txt');
%Error-free values
x = data(:,2); %X(time)values

%Vector of observations
L = data(:,3); %y(amplitude)values

%Number of observations
no_n = length(data);

%Initial values for the unknowns
% y = a+bx+cx^2-sin(x)
a=8.1977;
b=1.93665;
c=-0.39295;
 

%Vector of initial values for the unknowns
X_0 = [a b c]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n - no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
s_y = 0.05;        %[m]

%VC Matrix of the observations
S_LL = s_y^2 * eye(length(data));

%Theoretical standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_LL = sigma_0^-2*S_LL;

%Weight matrix
P = Q_LL^-1;
 
%--------------------------------------------------------------------------
%  Adjustment
% 
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-5;
delta = 10^-12;
max_x_hat =10^Inf; % starting value must be large 

%Number of iterations
iteration = 0;

while max_x_hat>epsilon || Check2>delta           
    
     %Observations as functions of the approximations for the unknowns
     L_0 = a+b*x+c*x.^2-sin(x);
     
     %Vector of reduced observations
     l = L-L_0;  
    
    % y = a+bx+cx^2-sin(x)
     %Design matrix with the elements from the Jacobian matrix J
     %      a             b       c               
     A =  [ones(6,1)      x     x.^2]; 
     
     %Normal matrix
     N = A'*P*A;
     
     %Vector of right hand side of normal equations
     n = A'*P*l;
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx = N^-1; 
    
     %Solution of the normal equations
     x_hat = Q_xx*n; 
       
     %Update
     X_0 = X_0+x_hat;
     a = X_0(1);
     b = X_0(2);
     c = X_0(3);
    
     %Check 1
     max_x_hat = max(abs(x_hat)); %--------------------------------------1
     
     %Vector of residuals
     v = A*x_hat-l; 
 
     %Vector of adjusted observations
     L_hat = L+v;
    
     %Objective function
     vTPv = v'*P*v; 
    
     %Functional relationships without the observations
     phi_X_hat = a+b*x+c*x.^2-sin(x);

     %Check 2
     Check2 = max(abs(L_hat-phi_X_hat)); %-------------------------------2
    
     %Update number of iterations
     iteration = iteration+1;
 end

if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

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
Q_vv = Q_LL - Q_LL_hat; 

%VC matrix of residuals
S_vv = s_0^2*Q_vv; 

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));





