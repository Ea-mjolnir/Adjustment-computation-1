%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%  Exercise 11: Adjustment Calculation - part VI  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 12, 2018
%   Last changes   : January 25, 2023
%
%--------------------------------------------------------------------------

clc;
close all;
format long g

%--------------------------------------------------------------------------
%   Task 1
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Coordinates - error free values
X_1=4527.150;      %[m]
Y_1=865.400 ;      %[m]
X_2=2047.250;      %[m]
Y_2=2432.550;      %[m]
X_3=27.150;        %[m]
Y_3=2865.22;       %[m]


S_100_1=6049.00;   %[m]
S_100_2=4736.830;  %[m]
S_100_3=5446.490;  %[m]
d_1_2 = sqrt((X_2-X_1)^2+(Y_2-Y_1)^2);                                     %[m]
b_1_2 = (atan2(Y_2-Y_1,X_2-X_1))*(180/3.142);                              %bearing in [degree]
alpha = acos((S_100_1^2+d_1_2^2-S_100_2^2)/(2*S_100_1*d_1_2))*(180/3.142); %included angle 
c =b_1_2-alpha; %bearing 1- point 100

%Vector of observations
L =[6049.00 4736.830 5446.490]'; %[m]

%Number of observations
no_n = length(L);

%Initial values for the unknowns
x100= 3728.9;                    %[m]
y100= 6861.5;                    %[m]

%Vector of initial values for the unknowns
X_0 =[x100 y100]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n-no_u;  

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
sig_1 = (1+(2*6.049))/1000;                     %[m]
sig_2 = (1+(2*4.736830))/1000;                  %[m]
sig_3 = (1+(2*5.446490))/1000;                  %[m]
S_LL =[sig_1^2;sig_2^2;sig_3^2].*eye(3);        %[m]

%Theoretical standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_LL = sigma_0^-2*S_LL;

%Weight matrix
P = Q_LL^-1;
%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-5;
delta = 10^-12;
max_x_hat =10^Inf;

%Number of iterations
iteration = 0;

while max_x_hat>epsilon || check2>delta          
    
     %Observations as functions of the approximations for the unknowns
     L_0(1,1)=sqrt((X_1-x100)^2 +(Y_1-y100)^2);
     L_0(2,1)=sqrt((X_2-x100)^2 +(Y_2-y100)^2);
     L_0(3,1)=sqrt((X_3-x100)^2 +(Y_3-y100)^2);
     %Vector of reduced observations
     l = L-L_0;
     
     %Design matrix with the elements from the Jacobian matrix J

      A(1,1)=-(X_1-x100)/(sqrt((X_1-x100)^2 +(Y_1-y100)^2));     
      A(1,2)=-(Y_1-y100)/(sqrt((X_1-x100)^2 +(Y_1-y100)^2));
      
      A(2,1)=-(X_2-x100)/(sqrt((X_2-x100)^2 +(Y_2-y100)^2));    
      A(2,2)=-(Y_2-y100)/(sqrt((X_2-x100)^2 +(Y_2-y100)^2));
      
      A(3,1)=-(X_3-x100)/(sqrt((X_3-x100)^2 +(Y_3-y100)^2));
      A(3,2)=-(Y_3-y100)/(sqrt((X_3-x100)^2 +(Y_3-y100)^2));
       
     %Normal matrix
       N = A'*P*A ; 
     
     %Vector of right hand side of normal equations
       n = A'*P*l;
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
       Q_xx = N^-1; 
    
     %Solution of the normal equations
       x_hat = Q_xx*n;
       
     %Update
       X_0 = X_0+x_hat; 
     
       x100 =  X_0(1);
       y100 =  X_0(2);
     
     %Check 1-----------------------------------------------------(1)
     max_x_hat = max(abs(x_hat));
          
     %Vector of residuals
     v = A*x_hat-l;
 
     %Vector of adjusted observations
     L_hat = L+v;
    
     %Objective function
     vTPv = v'*P*v;
    
     %Functional relationships without the observations
     phi_X_hat(1,1) = sqrt((X_1-x100)^2 +(Y_1-y100)^2);
     phi_X_hat(2,1) = sqrt((X_2-x100)^2 +(Y_2-y100)^2);
     phi_X_hat(3,1) = sqrt((X_3-x100)^2 +(Y_3-y100)^2);

     %Check 2---------------------------------------------------------------(2)
     check2 = max(abs(L_hat-phi_X_hat));
    
     %Update number of iterations
     iteration = iteration+1;
 
 end

if check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end


x100;
y100;

% Adjusted coordinate of point 100
Adj_100 = [x100 y100]';

%Empirical reference standard deviation
s_0 =sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx ;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));     

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat =s_0^2*Q_LL_hat; 

%Standard deviation of the adjusted observations
s_L_hat =sqrt(diag(S_LL_hat));

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));

disp('Observation => L  L_hat  s_L_hat  v  s_v')
obs = [L L_hat s_L_hat v s_v]
disp('Unknown => Adj_100 s_X')
unknown = [Adj_100 s_X]

