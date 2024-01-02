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
format short;

%--------------------------------------------------------------------------
%   Task 2: Adjustment of a parabola
%--------------------------------------------------------------------------  
%    y = ax^2+bx+c          ---------------- eqtn(parabola)
%   1.112  = a  + b  + c    ----------------  1
%   0.880  = 4a + 2b + c    ----------------  2
%   0.768  = 9a + 3b + c    ----------------  3
%   0.830  = 16a+ 4b + c    ----------------  4
%   1.175  = 25a+ 5b + c    ----------------  5
%--------------------------------------------------------------------------
%   Observations and redundancy
%--------------------------------------------------------------------------
%Vector of observations
L = [ 1.112 0.880 0.768 0.830 1.175]'; 
x = 1:5;
%Number of observations
no_n =5; 

%Number of unknowns
no_u =3; 

%Redundancy
r = no_n - no_u ;
%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
S_LL = (0.02^2.* eye(5));

%Theoretical reference standard deviation
sigma_0 = 0.01;         %a priori

%Cofactor matrix of the observations
Q_ll = ((sigma_0^2)^-1 * S_LL);
%Weight matrix
P = Q_ll^-1;
%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Design matrix
A =[1 1 1;
    4 2 1;
    9 3 1;
   16 4 1;
   25 5 1]; 
%Normal matrix
N = A'*P*A;        
%Vector of right hand side of normal equations
n = A'*P*L;    
%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_XX = N^-1;    
%Solution of normal equation
X_hat = Q_XX*n;
%Estimated unknown parameters
a = X_hat(1,1);
b = X_hat(2);
c = X_hat(3);     
%Vector of residuals
v =(A*X_hat)-L; 
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
s_0 = sqrt(vTPv/(no_n-no_u));      %a posteriori

%VC matrix of adjusted unknowns
E_XX_hat = (s_0^2*Q_XX);

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(E_XX_hat));

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_XX*A';

%VC matrix of adjusted observations
S_LL_hat = (s_0^2 * Q_LL_hat);
%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));
%Cofactor matrix of the residuals
Q_vv = Q_ll-Q_LL_hat;
%VC matrix of residuals
S_vv = s_0^2*Q_vv; 
%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));

%% ---------------------------task___5-------------------------------------------------
%Plot
L = [ 1.112 0.880 0.768 0.830 1.175]'; 
x = 1:5;
plot(x,L,'x') % 5 points
hold on
x_1= linspace(0,6);
plot(x_1,a*x_1.^2+ x_1.*b + c,'r')
hold off
xlabel('x-axis')
ylabel('y-axis')



