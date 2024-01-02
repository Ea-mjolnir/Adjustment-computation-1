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
%   Task 2
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Coordinates - error free values
X_1=4316.175;      %[m]
Y_1=935.411 ;      %[m]
X_2=4036.242;      %[m]
Y_2=2055.452;      %[m]
X_3=2136.262;      %[m]
Y_3=2331.535;      %[m]
X_4=1324.177;      %[m]
Y_4=1189.218;      %[m]
dist =[3491.901 2706.417 922.862 1819.298];  % distances [m]

d_1_4 =  sqrt((X_4-X_1)^2+(Y_4-Y_1)^2);                                    %[m]
b_1_4 = (atan2(Y_4-Y_1,X_4-X_1))*(180/3.142);                              %bearing in [degree]
alpha = acos((dist(1)^2+d_1_4^2-dist(4)^2)/(2*dist(1)*d_1_4))*(180/3.142); % included angle 

% Bearing 1-2
c = b_1_4 - alpha ;           % bearing of 1-P in degrees.
Px = X_1 +dist(1)*cosd(c);    % X initial coordinates of P  
Py = Y_1 +dist(1)*sind(c);    % Y initial coordinates of P
%Vector of observations
L =[3491.901 2706.417 922.862 1819.298]';     %[m]

%Number of observations
no_n = length(L);

%Initial values for the unknowns
P_x= 1500.67118774782;         %[m]
P_y= 3000.9159932146;          %[m]

%Vector of initial values for the unknowns
X_0 =[P_x P_y]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n-no_u;  

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
sig_1 = (2+(2*3.491901))/1000;  % in [m]
sig_2 = (2+(2*2.706417))/1000;  % in [m]
sig_3 = (2+(2*0.922862))/1000;  % in [m]
sig_4 = (2+(2*1.819298))/1000;  % in [m]

S_LL =[sig_1^2;sig_2^2;sig_3^2;sig_4].*eye(4);

%Theoretical standard deviation
sigma_0 = 0.01;     %a priori

%Cofactor matrix of the observations
Q_LL = sigma_0^-2*S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-5;
delta   = 10^-12;
max_x_hat =10^Inf;

%Number of iterations
iteration = 0;

while max_x_hat>epsilon || check_2>delta        
    
     %Observations as functions of the approximations for the unknowns
     
     L_0(1,1)= sqrt((X_1-P_x)^2 +(Y_1-P_y)^2);
     L_0(2,1)= sqrt((X_2-P_x)^2 +(Y_2-P_y)^2);
     L_0(3,1)= sqrt((X_3-P_x)^2 +(Y_3-P_y)^2);
     L_0(4,1)= sqrt((X_4-P_x)^2 +(Y_4-P_y)^2);
         
     %Vector of reduced observations
     l = L-L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
     
      A(1,1)= -(X_1-P_x)/(sqrt((X_1-P_x)^2 +(Y_1-P_y)^2));     
      A(1,2)= -(Y_1-P_y)/(sqrt((X_1-P_x)^2+(Y_1-P_y)^2));
      
      A(2,1)= -(X_2-P_x)/(sqrt((X_2-P_x)^2 +(Y_2-P_y)^2));     
      A(2,2)= -(Y_2-P_y)/(sqrt((X_2-P_x)^2 +(Y_2-P_y)^2));
      
      A(3,1)= -(X_3-P_x)/(sqrt((X_3-P_x)^2 +(Y_3-P_y)^2)); 
      A(3,2)= -(Y_3-P_y)/(sqrt((X_3-P_x)^2 +(Y_3-P_y)^2));
      
      A(4,1)= -(X_4-P_x)/(sqrt((X_4-P_x)^2 +(Y_4-P_y)^2));
      A(4,2)= -(Y_4-P_y)/(sqrt((X_4-P_x)^2 +(Y_4-P_y)^2));
      
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
     
      P_x =  X_0(1);
      P_y =  X_0(2);
     
     %Check 1-----------------------------------------------------------(1)
     max_x_hat = max(abs(x_hat));
          
     %Vector of residuals
     v = A*x_hat-l;
 
     %Vector of adjusted observations
     L_hat = L+v;
    
     %Objective function
     vTPv = v'*P*v;
    
     %Functional relationships without the observations
     phi_X_hat(1,1) = sqrt((X_1-P_x)^2 + (Y_1-P_y)^2);
     phi_X_hat(2,1) = sqrt((X_2-P_x)^2 + (Y_2-P_y)^2);
     phi_X_hat(3,1) = sqrt((X_3-P_x)^2 + (Y_3-P_y)^2);
     phi_X_hat(4,1) = sqrt((X_4-P_x)^2 + (Y_4-P_y)^2);

     %Check 2-----------------------------------------------------------(2)
      check_2 = max(abs(L_hat-phi_X_hat));
    
     %Update number of iterations
     iteration = iteration+1;
 
end
 
if check_2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

P_x; 
P_y; 

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
obs = [L  L_hat  s_L_hat  v  s_v]
disp('Unknown => Adj_100 s_X')
unknown = [Adj_100  s_X]
    

