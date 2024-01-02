%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%   Exercise 13: Adjustment Calculation - part VIII  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 15, 2018
%   Last changes   : February 08, 2023
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
% Observations and initial values for unknowns
%--------------------------------------------------------------------------
% Load all files
dist=load("Distances.txt");
dir=load("Directions.txt");
fixedpoint=load("FixedPoints.txt");
newpoint=load("NewPoints.txt");

% Vector of observations converted to radian
L =[dist(:,3);(dir(:,3)*pi/200)]; 

%Gauss-Krueger coordinates for control points [m]
for i=1:size(fixedpoint,1)
    eval(['y' num2str(fixedpoint(i,1)) '=' num2str(fixedpoint(i,2)) ';']);
    eval(['x' num2str(fixedpoint(i,1)) '=' num2str(fixedpoint(i,3)) ';']);
end

% Gauss-Krueger coordinates for initial values [m]
for i=1:size(newpoint,1)
    eval(['y' num2str(newpoint(i,1)) '=' num2str(newpoint(i,2)) ';']);
    eval(['x' num2str(newpoint(i,1)) '=' num2str(newpoint(i,3)) ';']);
end

% Initial values for orientation unknowns
w1 = 0;
w6 = 0;
w9 = 0;
w15 =0; 
% Initial values for unknowns
X_0 = [y1 x1 y15 x15 w1 w6 w9 w15]';

% Number of observations
no_n = length(L);

% Number of unknowns
no_u = length(X_0);
 
% Redundancy
r = no_n-no_u;


%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
% VC Matrix of the observations
sigma_a =0.1;            %[m] standard deviation for distances
sigma_b =0.001*pi/200;   %    standard deviation in radians
u =repelem(sigma_a^2,5);
v =repelem(sigma_b^2,9);

S_LL=[u,v]'.*eye(14);   %Vc matrix

% Theoretical standard deviation
sigma_0 = 1;

% Cofactor matrix of the observations
Q_LL = 1/sigma_0^2*S_LL;

% Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
% break-off condition
epsilon = 10^-5; 
delta = 10^-9;
max_x_hat = 10^Inf;

% Number of iterations
iteration = 0;

% Initialising A
A = zeros(no_n,no_u);

% Iteration
while max_x_hat>epsilon || Check2>delta

	% Distances
	L_0(1) =dis(y6,x6,y1,x1);
    L_0(2) =dis(y9,x9,y1,x1); 
	L_0(3) =dis(y9,x9,y6,x6);
    L_0(4) =dis(y15,x15,y1,x1);
    L_0(5) =dis(y15,x15,y9,x9);
    
	% Directions
	L_0(6)=  direction(y1,x1,y6,x6,w1);
    L_0(7)=  direction(y1,x1,y15,x15,w1);
    L_0(8)=  direction(y6,x6,y1,x1,w6);
    L_0(9)=  direction(y6,x6,y9,x9,w6);
    L_0(10)= direction(y9,x9,y15,x15,w9);
    L_0(11)= direction(y9,x9,y1,x1,w9);
    L_0(12)= direction(y9,x9,y6,x6,w9);
    L_0(13)= direction(y15,x15,y1,x1,w15);
    L_0(14)= direction(y15,x15,y9,x9,w15);

    % Vector of reduced observations
    l = L-L_0';

    % Design matrix with the elements from the Jacobian matrix J
    % y1 x1 y15 x15 w1 w6 w9 w15
	A(1,1) = ds_dy_to(y6,x6,y1,x1);
    A(1,2) = ds_dx_to(y6,x6,y1,x1);
    
    A(2,1) = ds_dy_to(y9,x9,y1,x1);
    A(2,2) = ds_dx_to(y9,x9,y1,x1);
    
    A(4,1) = ds_dy_to(y15,x15,y1,x1);
    A(4,2) = ds_dx_to(y15,x15,y1,x1);
    A(4,3) = ds_dy_from(y15,x15,y1,x1);
    A(4,4) = ds_dx_from(y15,x15,y1,x1);
    
    A(5,3) = ds_dy_from(y15,x15,y9,x9);
    A(5,4) = ds_dx_from(y15,x15,y9,x9);
    
    A(6,1) = dr_dy_from(y1,x1,y6,x6);
    A(6,2) = dr_dx_from(y1,x1,y6,x6);
    A(6,5) = -1;
    
    A(7,1) = dr_dy_from(y1,x1,y15,x15);
    A(7,2) = dr_dx_from(y1,x1,y15,x15);
    A(7,3) = dr_dy_to(y1,x1,y15,x15);
    A(7,4) = dr_dx_to(y1,x1,y15,x15);
    A(7,5) = -1;
    
    A(8,1) = dr_dy_to(y6,x6,y1,x1);
    A(8,2) = dr_dx_to(y6,x6,y1,x1);
    A(8,6) = -1;
    
    A(9,6) = -1;
    
    A(10,3) = dr_dy_to(y9,x9,y15,x15);
    A(10,4) = dr_dx_to(y9,x9,y15,x15);
    A(10,7) = -1;
    
    A(11,1) = dr_dy_to(y9,x9,y1,x1);
    A(11,2) = dr_dx_to(y9,x9,y1,x1);
    A(11,7) = -1;
    
    A(12,7) = -1;
    
    A(13,1) = dr_dy_to(y15,x15,y1,x1);
    A(13,2) = dr_dx_to(y15,x15,y1,x1);
    A(13,3) = dr_dy_from(y15,x15,y1,x1);
    A(13,4) = dr_dx_from(y15,x15,y1,x1);
    A(13,8) = -1;
    
    A(14,3) = dr_dy_from(y15,x15,y9,x9);
    A(14,4) = dr_dx_from(y15,x15,y9,x9);
    A(14,8) = -1;
    
    % Normal matrix
    N = A'*P*A;

    % Vector of right hand side of normal equations
    n = A'*P*l;

    % Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx = inv(N);
    
    % Solution of normal equations
    x_hat = Q_xx*n;
    
    % Adjusted unknowns
    X_hat = X_0+x_hat;
    
    % Update
    X_0 = X_hat;

    y1 =X_0(1);
    x1 =X_0(2);
    y15=X_0(3);
    x15=X_0(4);
    w1 =X_0(5);
    w6 =X_0(6);
    w9 =X_0(7);
    w15=X_0(8);
    
    % Check 1
    max_x_hat = max(abs(x_hat));
  
    % Vector of residuals
    v = A*x_hat-l;
 
    % Vector of adjusted observations
    L_hat= L+v;
    
    % Function
    vTPv = v'*P*v;
    
    % Functional relationships 
    phi_X_hat  =[dis(y6,x6,y1,x1)
                 dis(y9,x9,y1,x1) 
	             dis(y9,x9,y6,x6)
                 dis(y15,x15,y1,x1);
                 dis(y15,x15,y9,x9);
                 direction(y1,x1,y6,x6,w1)
                 direction(y1,x1,y15,x15,w1)
                 direction(y6,x6,y1,x1,w6)
                 direction(y6,x6,y9,x9,w6)
                 direction(y9,x9,y15,x15,w9)
                 direction(y9,x9,y1,x1,w9)
                 direction(y9,x9,y6,x6,w9)
                 direction(y15,x15,y1,x1,w15)
                 direction(y15,x15,y9,x9,w15)];
	 
    %Check 2
    Check2 = max(abs(L_hat-phi_X_hat)); 
    
    % Update number of iterations
    iteration = iteration+1;

end

if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

% Convert to [gon]
for i=1:length(X_0)
    while i>4
        if X_0(i)<0
            X_0(i)=(X_0(i)*200/pi)+400;
        else
            X_0(i)=X_0(i)*200/pi;
        end
     break
    end    
end

X= ["y1:" "x1:" "y15:" "x15:" "w1:" "w6:" "w9:" "w15:"]';
Results=[X,X_0];

% Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

% VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

% Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));
for i=1:length(s_X)
    while i>4
        s_X(i)=s_X(i)*200/pi;
       break
    end
end
std_unknowns =[Results,s_X];

% Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

% VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

% Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));
max_len = max(length(s_L_hat),length(L_hat));
for i = 1:max_len
    if i > 5 % changed while loop to if statement, since the loop should only execute once
        s_L_hat(i) = s_L_hat(i)*200/pi;
        L_hat(i) = L_hat(i)*(200/pi);
        break
    end
end

std_obs =[L,L_hat,s_L_hat];
        
% Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

% VC matrix of residuals
S_vv = s_0^2*Q_vv;              

% Standard deviation of the residuals
s_v = sqrt(diag(S_vv));

for i=1:length(s_v)
    while i>5
        s_v(i)=s_v(i)*200/pi;
        break
    end
end
std_resid =[v,s_v];



