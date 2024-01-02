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
close all;
format long g;

%--------------------------------------------------------------------------
%   Task 1
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Load files
 dir  =load("directions.txt");   %[gon]
coord =load("Points.txt");     %[m]   %Error-free

for i=1:size(coord,1)
    eval(['y' num2str(coord(i,1)) '=' num2str(coord(i,2)) ';']);
    eval(['x' num2str(coord(i,1)) '=' num2str(coord(i,3)) ';']);
end

%Vector of observations
L =dir(:,3)*(pi/200); % in rad

%Number of observations
no_n = length(L);

%Initial values for the unknowns
 x3 =343.623;   %[m]
 y3 =236.622;   %[m]
 w3 =0;

%Vector of initial values for the unknowns
X_0 =[x3 y3 w3]'; 

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n-no_u;  

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
s_dir = 0.001*(pi/200);   %[gon]->[rad]

%VC Matrix of the observations
S_LL = (s_dir^2)*eye(5); 

%Theoretical standard deviation
sigma_0 = 1;     %a priori

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
max_x_hat =10^Inf; 

%Number of iterations
iteration = 0;

while max_x_hat>epsilon || check2>delta           
    
     %Observations as functions of the approximations for the unknowns
     L_0(1)=atan2(y3-y1,x3-x1);
     L_0(2)=atan2(y3-y2,x3-x2);
     L_0(3)=atan2(y3-y4,x3-x4);
     L_0(4)=atan2(y3-y5,x3-x5);
     L_0(5)=atan2(y3-y6,x3-x6);
     
     for i =1:length(L_0)
         if L_0(i)<0
             L_0(i) =(L_0(i)+(2*pi))-w3;
         else
             L_0(i)=L_0(i)-w3;
         end
     end
     %Vector of reduced observations
     l =L-L_0';
    
     %Design matrix with the elements from the Jacobian matrix J
     A(1,1)= (y1-y3)/(dis(y3,x3,y1,x1))^2;
     A(1,2)=-(x1-x3)/(dis(y3,x3,y1,x1))^2;
     A(1,3)=-1;

     A(2,1)= (y2-y3)/(dis(y3,x3,y2,x2))^2;
     A(2,2)=-(x2-x3)/(dis(y3,x3,y2,x2))^2;
     A(2,3)=-1;

     A(3,1)= (y4-y3)/(dis(y3,x3,y4,x4))^2;
     A(3,2)=-(x4-x3)/(dis(y3,x3,y4,x4))^2;
     A(3,3)=-1;

     A(4,1)= (y5-y3)/(dis(y3,x3,y5,x5))^2;
     A(4,2)=-(x5-x3)/(dis(y3,x3,y5,x5))^2;
     A(4,3)=-1;

     A(5,1)= (y6-y3)/(dis(y3,x3,y6,x6))^2;
     A(5,2)=-(x6-x3)/(dis(y3,x3,y6,x6))^2;
     A(5,3)=-1;

     
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
     w3 = X_0(3);

    
     %Check 1-----------------------------------------------(1)
     max_x_hat = max(abs(x_hat));
     
     %Vector of residuals
     v =A*x_hat - l; 
 
     %Vector of adjusted observations
     L_hat = L+v;
    
     %Function
     vTPv = v'*P*v;
    
     %Functional relationships without the observations
     phi_X_hat= [atan2(y3-y1,x3-x1)
                 atan2(y3-y2,x3-x2)
                 atan2(y3-y4,x3-x4)
                 atan2(y3-y5,x3-x5)
                 atan2(y3-y6,x3-x6)];
     for i =1:length(phi_X_hat)
         if phi_X_hat(i)<0
             phi_X_hat(i) =(phi_X_hat(i)+(2*pi))-w3;
         else
             phi_X_hat(i)=phi_X_hat(i)-w3;
         end
     end

     %Check 2-------------------------------------------------(2)
     check2 = max(abs(L_hat-phi_X_hat));
    
     %Update number of iterations
     iteration = iteration+1;
  
end

if check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

x3
y3
w3

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));        %[m]

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

    
