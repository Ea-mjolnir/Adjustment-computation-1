%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%    Exercise 12: Adjustment Calculation - part VII  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 12, 2018
%   Last changes   : January 31, 2023
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
%Load files
dir   =load("directions.txt");   %[gon]
coord =load("Points.txt");       %[m]   %Error-free

for i=1:size(coord,1)
   eval(['y' num2str(coord(i,1)) '=' num2str(coord(i,2)) ';']);
   eval(['x' num2str(coord(i,1)) '=' num2str(coord(i,3)) ';']);
end


%Vector of observations
L =dir(:,3)*pi/200; % in rad
L_Ang =[L(3)-L(2);L(4)-L(3);L(5)-L(4);L(1)-L(5)]; % included angles in [rad]
%Number of observations
no_n = length(L_Ang);
 
 % initial values of point x3 
 x3 =343.623; %[m]
 y3 =236.622; %[m]

%Vector of initial values for the unknowns
X_0 =[x3 y3]'; 

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n-no_u;  

% --------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
s_dir = 0.001*(pi/200); %[gon]->[rad]

%VC Matrix of the observations
F =[ 0 -1  1  0  0
     0  0 -1  1  0
     0  0  0 -1  1
     1  0  0  0 -1];

S_LL = (s_dir^2)*eye(5);
S_xx =  F*S_LL*F';

%Theoretical standard deviation
sigma_0 = 1;     %a priori

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_xx;

%Weight matrix
P = inv(Q_LL);

% ------------------------------------------------------------------------
%    Adjustment
%  ------------------------------------------------------------------------
%break-off conditions
epsilon= 10^-5;
delta = 10^-12;
max_x_hat=10^Inf; 

%Number of iterations
iteration = 0;

while all(max_x_hat>epsilon) || all(Check2>delta)            
    
     %Observations as functions of the approximations for the unknowns in
     %radians
     L_0(1)=atan2(y4-y3,x4-x3)-atan2(y2-y3,x2-x3);
     L_0(2)=atan2(y5-y3,x5-x3)-atan2(y4-y3,x4-x3);
     L_0(3)=atan2(y6-y3,x6-x3)-atan2(y5-y3,x5-x3);
     L_0(4)=atan2(y1-y3,x1-x3)-atan2(y6-y3,x6-x3);
     
     for i =1:length(L_0)
         if L_0(i)<0
             L_0(i)=(L_0(i)+(2*pi));
         else
             L_0(i)=L_0(i);
         end
     end
     %Vector of reduced observations
     l =L_Ang - L_0';
    
     %Design matrix with the elements from the Jacobian matrix J
     A(1,1)=  (y4-y3)/dis(x3,y3,x4,y4)^2 - (y2-y3)/dis(x3,y3,x2,y2)^2;
     A(1,2)= -(x4-x3)/dis(x3,y3,x4,y4)^2 + (x2-x3)/dis(x3,y3,x2,y2)^2;
  

     A(2,1)=  (y5-y3)/dis(x3,y3,x5,y5)^2 - (y4-y3)/dis(x3,y3,x4,y4)^2;
     A(2,2)= -(x5-x3)/dis(x3,y3,x5,y5)^2 + (x4-x3)/dis(x3,y3,x4,y4)^2;
  

     A(3,1)=  (y6-y3)/dis(x3,y3,x6,y6)^2 - (y5-y3)/dis(x3,y3,x5,y5)^2;
     A(3,2)= -(x6-x3)/dis(x3,y3,x6,y6)^2 + (x5-x3)/dis(x3,y3,x5,y5)^2;
    

     A(4,1)=  (y1-y3)/dis(x3,y3,x1,y1)^2 - (y6-y3)/dis(x3,y3,x6,y6)^2;
     A(4,2)= -(x1-x3)/dis(x3,y3,x1,y1)^2 + (x6-x3)/dis(x3,y3,x6,y6)^2;
    
     
     %Normal matrix
     N =A'*P*A;
     
     %Vector of right hand side of normal equations
     n = A'*P*l;
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx = N^-1;
    
     %Solution of the normal equations
     x_hat = Q_xx*n;
       
     %Update
     X_0 = X_0+x_hat;
   
     x3 = X_0(1);
     y3 = X_0(2);
    
     
     %Check 1-----------------------------------------------(1)
     max_x_hat = max(abs(x_hat));
     
     %Vector of residuals
     v =A*x_hat - l; 
 
     %Vector of adjusted observations
     L_hat = L_Ang+v;
    
     %Function
     vTPv = v'*P*v;
    
     %Functional relationships without the observations
     
     phi_X_hat=[atan2(y4-y3,x4-x3)-atan2(y2-y3,x2-x3)
                atan2(y5-y3,x5-x3)-atan2(y4-y3,x4-x3)
                atan2(y6-y3,x6-x3)-atan2(y5-y3,x5-x3)
                atan2(y1-y3,x1-x3)-atan2(y6-y3,x6-x3)];
                  
     for i =1:length(phi_X_hat)
         if phi_X_hat(i)<0
             phi_X_hat(i) =(phi_X_hat(i)+(2*pi));
         else
             phi_X_hat(i)=phi_X_hat(i);
         end
     end
     
     %Check 2-------------------------------------------------(2)
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


%Convert to [gon]
L_Ang_gon = L_Ang *(200/pi);
L_hat_gon = L_hat*200/pi;
v_gon = v*(200/pi);

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat)); %[m]

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));     
s_L_hat_gon =s_L_hat*(200/pi); % back to gon

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv)); 
s_v_gon =s_v*(200/pi);

disp('Observation =>L_Ang_gon ,L_hat_gon,s_L_hat_gon,v_gon,s_v_gon')
obs = [L_Ang_gon,L_hat_gon,s_L_hat_gon,v_gon,s_v_gon]
disp('Unknown => X3 Y3 & s_X ')
unknown = [X_0,s_X] 
