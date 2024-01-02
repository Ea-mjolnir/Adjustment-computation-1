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
% Load files
dir = load("directions.txt"); %[gon]
coord =load("Points7.txt");   %[m]     %Error-free

%Vector of observations
L = dir(:,3)*pi/200;    %[gon]->[rad]
L_gon = dir(:,3); %[gon]

for i=1:size(coord,1)
    eval(['y' num2str(coord(i,1)) '=' num2str(coord(i,2)) ';']);
    eval(['x' num2str(coord(i,1)) '=' num2str(coord(i,3)) ';']);
end

%Number of observations
no_n = length(L);

%Initial values for the unknowns
x3 = 242.900; 
y3 = 493.700;
w3 = 0;

%Vector of initial values for the unknowns
X_0 = [x3 y3 w3]';

%Number of unknowns
no_u = length(X_0);

%Number of constraints
no_b = 1;

%Redundancy
r = no_n - no_u + no_b;   

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
s_dir = 0.001 * pi/200;   %[gon]->[rad]

%VC Matrix of the observations
S_LL = s_dir^2*eye(5);   

%Theoretical standard deviation
sigma_0 = 1;  %a priori

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-5;
delta = 10^-12;
max_x_hat = 10^Inf;

%Number of iterations
iteration = 0;
A=zeros(5,3);
while all(max_x_hat>epsilon) ||all(Check2>delta)           
    
     %Observations as functions of the approximations for the unknowns
     L_0(1) = direction(y3,x3,y1,x1,w3);
     L_0(2) = direction(y3,x3,y2,x2,w3);
     L_0(3) = direction(y3,x3,y4,x4,w3);
     L_0(4) = direction(y3,x3,y5,x5,w3);
     L_0(5) = direction(y3,x3,y6,x6,w3);
    %L_0(6) = direction(y3,x3,y7,x7,w3);
     %Vector of reduced observations
     l = L-L_0';
    
     %Design matrix A with the elements from the Jacobian matrix J
        %       x3    y3    w3
      A(1,1) = dr_dx_from(y3,x3,y1,x1);   % bearing 3 --> 1
      A(1,2) = dr_dy_from(y3,x3,y1,x1);
      A(1,3) = -1;

      A(2,1) = dr_dx_from(y3,x3,y2,x2);   % bearing 3 --> 2
      A(2,2) = dr_dy_from(y3,x3,y2,x2);
      A(2,3) = -1;

      A(3,1) = dr_dx_from(y3,x3,y4,x4);   % bearing 3 --> 4
      A(3,2) = dr_dy_from(y3,x3,y4,x4);
      A(3,3) = -1;      
      
      A(4,1) = dr_dx_from(y3,x3,y5,x5);   % bearing 3 --> 5
      A(4,2) = dr_dy_from(y3,x3,y5,x5);
      A(4,3) = -1;
     
      A(5,1) = dr_dx_from(y3,x3,y6,x6);   % bearing 3 --> 6
      A(5,2) = dr_dy_from(y3,x3,y6,x6);
      A(5,3) = -1;

      %A(6,1) = dr_dx_from(y3,x3,y7,x7);
      %A(6,2) = dr_dy_from(y3,x3,y7,x7);
      %A(6,3) = -1;
     
     %Design matrix C with the elements from the Jacobian matrix J
     %                         x         
     C = [-2*(x7-X_0(1))/(2*((x7-X_0(1))^2+(y7-X_0(2))^2)^(1/2))   
          -2*(y7-X_0(2))/(2*((x7-X_0(1))^2+(y7-X_0(2))^2)^(1/2))    
                               0                                ];%constraint

     %Normal matrix
     N = A'*P*A ;
     
     %Extended normal matrix
     N_ext = [N   C 
              C'  0 ];
     
     %Vector of right hand side of normal equations
     n = A'*P*l;
     
     %Extended vector of right hand side of normal equations
     
     n_ext = [                 n  
              25-sqrt((x7-X_0(1))^2+(y7-X_0(2))^2)];   %constraint
          
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx_ext = N_ext^-1;
     Q_xx = Q_xx_ext(1:no_u,1:no_u);

     %Solution of the normal equations
     x_hat = Q_xx_ext*n_ext;
       
     %Update
     X_0 = X_0 + x_hat(1:no_u);
    
     x3 = X_0(1);
     y3 = X_0(2);
     w3 = X_0(3);

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
      phi_X_hat = [direction(y3,x3,y1,x1,w3);
                   direction(y3,x3,y2,x2,w3);
                   direction(y3,x3,y4,x4,w3);
                   direction(y3,x3,y5,x5,w3);
                   direction(y3,x3,y6,x6,w3)];
                  

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


x3
y3
w3=(w3*200/pi)+400
s=sqrt((x7-x3)^2+(y7-y3)^2)

%Convert to [gon]
L_hat_gon = L_hat*200/pi;
v_gon = v*200/pi;
v_gon(1) = v_gon(1)+400;
v_gon(2) = v_gon(2)+400;
X_0(3) = (X_0(3)*200/pi) +400;

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
       
s_X = sqrt(diag(S_XX_hat));        %[m]
s_X_gon = s_X(3)*200/pi;           %[gon] for s_w3


%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));     
s_L_hat_gon = s_L_hat*200/pi;      %[rad]->[gon]

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));         
s_v_gon = s_v*200/pi;              %[rad]->[gon]

disp('Observation => L  L_hat  s_L_hat  v  s_v ')
obs = [L_gon,L_hat_gon,s_L_hat_gon,v_gon,s_v_gon]
disp('Unknown => X3 Y3 w3 & s_X ')
unknown = [X_0 , s_X]         
    

