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
Cp=load("GK control points.txt");
initial_Cp=load("GK_initial _controlpoints.txt");
 
% Vector of observations converted to radian
L =[dist(:,3);(dir(:,3)*pi/200)]; 

%Gauss-Krueger coordinates for control points [m]
for i=1:size(Cp,1)
    eval(['y' num2str(Cp(i,1)) '=' num2str(Cp(i,2)) ';']);
    eval(['x' num2str(Cp(i,1)) '=' num2str(Cp(i,3)) ';']);
end

% Gauss-Krueger coordinates for initial values [m]
for i=1:size(initial_Cp,1)
    eval(['y' num2str(initial_Cp(i,1)) '=' num2str(initial_Cp(i,2)) ';']);
    eval(['x' num2str(initial_Cp(i,1)) '=' num2str(initial_Cp(i,3)) ';']);
end
 
% Initial values for orientation unknowns
w1000 =0;
w2000 =0;
w3000 =0;
w100 =0; 
w101 =0; 
w102 =0;
w103 =0;
% Initial values for unknowns
X_0 = [y100 x100 y101 x101 y102 x102 y103 x103 w1000 w2000 w3000 w100 w101 w102 w103]';

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
sigma_a = 0.001;           %[m] standard deviation for distances
sigma_b =0.001*pi/200;     % standard deviation in radians
u =repelem(sigma_a^2,18);
v =repelem(sigma_b^2,18);

S_LL =[u,v]'.*eye(36);  %Vc matrix

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
	L_0(1) =dis(y1000,x1000,y100,x100);
    L_0(2) =dis(y1000,x1000,y102,x102); 
	L_0(3) =dis(y2000,x2000,y103,x103);
    L_0(4) =dis(y2000,x2000,y101,x101);
    L_0(5) =dis(y3000,x3000,y100,x100);
    L_0(6) =dis(y3000,x3000,y103,x103);
    L_0(7) =dis(y100,x100,y1000,x1000);
    L_0(8) =dis(y100,x100,y102,x102);
    L_0(9) =dis(y100,x100,y2000,x2000);
    L_0(10)=dis(y100,x100,y3000,x3000);
    L_0(11)=dis(y101,x101,y103,x103);
    L_0(12)=dis(y101,x101,y2000,x2000);
    L_0(13)=dis(y102,x102,y1000,x1000);
    L_0(14)=dis(y102,x102,y100,x100);
    L_0(15)=dis(y102,x102,y2000,x2000);
    L_0(16)=dis(y103,x103,y3000,x3000);
    L_0(17)=dis(y103,x103,y2000,x2000);
    L_0(18)=dis(y103,x103,y101,x101);
    
	% Directions
	L_0(19)=direction(y1000,x1000,y100,x100,w1000);
    L_0(20)=direction(y1000,x1000,y102,x102,w1000); 
	L_0(21)=direction(y2000,x2000,y103,x103,w2000);
    L_0(22)=direction(y2000,x2000,y101,x101,w2000);
    L_0(23)=direction(y3000,x3000,y100,x100,w3000);
    L_0(24)=direction(y3000,x3000,y103,x103,w3000);
    L_0(25)=direction(y100,x100,y1000,x1000,w100);
    L_0(26)=direction(y100,x100,y102,x102,w100);
    L_0(27)=direction(y100,x100,y2000,x2000,w100);
    L_0(28)=direction(y100,x100,y3000,x3000,w100);
    L_0(29)=direction(y101,x101,y103,x103,w101);
    L_0(30)=direction(y101,x101,y2000,x2000,w101);
    L_0(31)=direction(y102,x102,y1000,x1000,w102);
    L_0(32)=direction(y102,x102,y100,x100,w102);
    L_0(33)=direction(y102,x102,y2000,x2000,w102);
    L_0(34)=direction(y103,x103,y3000,x3000,w103);
    L_0(35)=direction(y103,x103,y2000,x2000,w103);
    L_0(36)=direction(y103,x103,y101,x101,w103);
    % Vector of reduced observations
    l = L-L_0';

% Derivative of the observations with respect to the unknowns 

% Design matrix with the elements from the Jacobian matrix J
% 1:y100 2:x100 3:y101 4:x101 5:y102 6:x102 7:y103 8:x103 9:w1000 10:w2000 11:w3000 12:w100 13:w101 14:w102 15:w103
	A(1,1) = ds_dy_to(y1000,x1000,y100,x100);
    A(1,2) = ds_dx_to(y1000,x1000,y100,x100);
    
    A(2,5) = ds_dy_to(y1000,x1000,y102,x102);
    A(2,6) = ds_dx_to(y1000,x1000,y102,x102);
    
    A(3,7) = ds_dy_to(y2000,x2000,y103,x103);
    A(3,8) = ds_dx_to(y2000,x2000,y103,x103);
    
    A(4,3) = ds_dy_to(y2000,x2000,y101,x101);
    A(4,4) = ds_dx_to(y2000,x2000,y101,x101);
    
    A(5,1) = ds_dy_to(y3000,x3000,y100,x100);
    A(5,2) = ds_dx_to(y3000,x3000,y100,x100);
    
    A(6,7) = ds_dy_to(y3000,x3000,y103,x103);
    A(6,8) = ds_dx_to(y3000,x3000,y103,x103);
    
    A(7,1) = ds_dy_from(y100,x100,y1000,x1000);
    A(7,2) = ds_dx_from(y100,x100,y1000,x1000);
    
%1:y100 2:x100 3:y101 4:x101 5:y102 6:x102 7:y103 8:x103 9:w1000 10:w2000 11:w3000 12:w100 13:w101 14:w102 15:w103    
    A(8,1) = ds_dy_from(y100,x100,y102,x102);
    A(8,2) = ds_dx_from(y100,x100,y102,x102);
    A(8,5) = ds_dy_to(y100,x100,y102,x102);
    A(8,6) = ds_dx_to(y100,x100,y102,x102); 
    
    A(9,1) = ds_dy_from(y100,x100,y2000,x2000);
    A(9,2) = ds_dx_from(y100,x100,y2000,x2000);
    
    A(10,1) = ds_dy_from(y100,x100,y3000,x3000);
    A(10,2) = ds_dx_from(y100,x100,y3000,x3000); 
%
    A(11,3) = ds_dy_from(y101,x101,y103,x103);
    A(11,4) = ds_dx_from(y101,x101,y103,x103); 
    A(11,7) = ds_dy_to(y101,x101,y103,x103);
    A(11,8) = ds_dx_to(y101,x101,y103,x103);
    
    A(12,3) = ds_dy_from(y101,x101,y2000,x2000);
    A(12,4) = ds_dx_from(y101,x101,y2000,x2000);    
%1:y100 2:x100 3:y101 4:x101 5:y102 6:x102 7:y103 8:x103 9:w1000 10:w2000 11:w3000 12:w100 13:w101 14:w102 15:w103      
    A(13,5) = ds_dy_from(y102,x102,y1000,x1000);
    A(13,6) = ds_dx_from(y102,x102,y1000,x1000);
    
    A(14,5) = ds_dy_from(y102,x102,y100,x100);
    A(14,6) = ds_dx_from(y102,x102,y100,x100);
    A(14,1) = ds_dy_to(y102,x102,y100,x100);
    A(14,2) = ds_dx_to(y102,x102,y100,x100);
    
    A(15,5) = ds_dy_from(y102,x102,y2000,x2000);
    A(15,6) = ds_dx_from(y102,x102,y2000,x2000); 
    
    A(16,7) = ds_dy_from(y103,x103,y3000,x3000);
    A(16,8) = ds_dx_from(y103,x103,y3000,x3000); 
    
    A(17,7) = ds_dy_from(y103,x103,y2000,x2000);
    A(17,8) = ds_dx_from(y103,x103,y2000,x2000);
    
    A(18,7) = ds_dy_from(y103,x103,y101,x101);
    A(18,8) = ds_dx_from(y103,x103,y101,x101); 
    A(18,3) = ds_dy_to(y103,x103,y101,x101);
    A(18,4) = ds_dx_to(y103,x103,y101,x101);
 
%1:y100 2:x100 3:y101 4:x101 5:y102 6:x102 7:y103 8:x103 9:w1000 10:w2000 11:w3000 12:w100 13:w101 14:w102 15:w103      

 % For directions 
    A(19,1) = dr_dy_to(y1000,x1000,y100,x100);
    A(19,2) = dr_dx_to(y1000,x1000,y100,x100);
    A(19,9) = -1;
    
    A(20,5) = dr_dy_to(y1000,x1000,y102,x102);
    A(20,6) = dr_dx_to(y1000,x1000,y102,x102);
    A(20,9) = -1;
    
    A(21,7) = dr_dy_to(y2000,x2000,y103,x103);
    A(21,8) = dr_dx_to(y2000,x2000,y103,x103);
    A(21,10) = -1;
%1:y100 2:x100 3:y101 4:x101 5:y102 6:x102 7:y103 8:x103 9:w1000 10:w2000 11:w3000 12:w100 13:w101 14:w102 15:w103     
    A(22,3) = dr_dy_to(y2000,x2000,y101,x101);
    A(22,4) = dr_dx_to(y2000,x2000,y101,x101);
    A(22,10) = -1;    
    
    A(23,1) = dr_dy_to(y3000,x3000,y100,x100);
    A(23,2) = dr_dx_to(y3000,x3000,y100,x100);
    A(23,11) = -1;  
    
    A(24,7) = dr_dy_to(y3000,x3000,y103,x103);
    A(24,8) = dr_dx_to(y3000,x3000,y103,x103);
    A(24,11) = -1;
    
    A(25,1) = dr_dy_from(y100,x100,y1000,x1000);
    A(25,2) = dr_dx_from(y100,x100,y1000,x1000);
    A(25,12) = -1;    
%1:y100 2:x100 3:y101 4:x101 5:y102 6:x102 7:y103 8:x103 9:w1000 10:w2000 11:w3000 12:w100 13:w101 14:w102 15:w103         
    A(26,1) = dr_dy_from(y100,x100,y102,x102);
    A(26,2) = dr_dx_from(y100,x100,y102,x102);
    A(26,5) = dr_dy_to(y100,x100,y102,x102);
    A(26,6) = dr_dx_to(y100,x100,y102,x102);
    A(26,12)= -1;
    A(27,1) = dr_dy_from(y100,x100,y2000,x2000);
    A(27,2) = dr_dx_from(y100,x100,y2000,x2000);
    A(27,12)= -1;   
    A(28,1) = dr_dy_from(y100,x100,y3000,x3000);
    A(28,2) = dr_dx_from(y100,x100,y3000,x3000);
    A(28,12)= -1;    
    A(29,3) = dr_dy_from(y101,x101,y103,x103);
    A(29,4) = dr_dx_from(y101,x101,y103,x103);
    A(29,7) = dr_dy_to(y101,x101,y103,x103);
    A(29,8) = dr_dx_to(y101,x101,y103,x103);
    A(29,13)= -1;
%1:y100 2:x100 3:y101 4:x101 5:y102 6:x102 7:y103 8:x103 9:w1000 10:w2000 11:w3000 12:w100 13:w101 14:w102 15:w103        
    A(30,3) = dr_dy_from(y101,x101,y2000,x2000);
    A(30,4) = dr_dx_from(y101,x101,y2000,x2000);
    A(30,13)= -1;
    A(31,5) = dr_dy_from(y102,x102,y1000,x1000);
    A(31,6) = dr_dx_from(y102,x102,y1000,x1000);
    A(31,14)= -1;
    
    A(32,5) = dr_dy_from(y102,x102,y100,x100);
    A(32,6) = dr_dx_from(y102,x102,y100,x100);
    A(32,1) = dr_dy_to(y102,x102,y100,x100);
    A(32,2) = dr_dx_to(y102,x102,y100,x100);   
    A(32,14)= -1;
%1:y100 2:x100 3:y101 4:x101 5:y102 6:x102 7:y103 8:x103 9:w1000 10:w2000 11:w3000 12:w100 13:w101 14:w102 15:w103       
    A(33,5) = dr_dy_from(y102,x102,y2000,x2000);
    A(33,6) = dr_dx_from(y102,x102,y2000,x2000);
    A(33,14) = -1;
    
    A(34,7) = dr_dy_from(y103,x103,y3000,x3000);
    A(34,8) = dr_dx_from(y103,x103,y3000,x3000);
    A(34,15) = -1;    

    A(35,7) = dr_dy_from(y103,x103,y2000,x2000);
    A(35,8) = dr_dx_from(y103,x103,y2000,x2000);
    A(35,15) = -1; 
    
    A(36,7) = dr_dy_from(y103,x103,y101,x101);
    A(36,8) = dr_dx_from(y103,x103,y101,x101);
    A(36,3) = dr_dy_to(y103,x103,y101,x101);
    A(36,4) = dr_dx_to(y103,x103,y101,x101);  
    A(36,15)= -1; 
    
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
%X_0 = [y100 x100 y101 x101 y102 x102 y103 x103 w1000 w2000 w3000 w100 w101 w102 w103]';
    y100 =X_0(1);
    x100 =X_0(2);
    y101 =X_0(3);
    x101 =X_0(4);  
    y102 =X_0(5);
    x102 =X_0(6);    
    y103 =X_0(7);
    x103 =X_0(8);
    w1000 =X_0(9);
    w2000 =X_0(10);
    w3000 =X_0(11);
    w100 =X_0(12);
    w101 =X_0(13);
    w102 =X_0(14);
    w103 =X_0(15);
    
    % Check 1
    max_x_hat = max(abs(x_hat));
  
    % Vector of residuals
    v = A*x_hat-l;
 
    % Vector of adjusted observations
    L_hat= L+v;
    
    % Function
    vTPv = v'*P*v;
    
    % Functional relationships 
    phi_X_hat =[dis(y1000,x1000,y100,x100)
                dis(y1000,x1000,y102,x102) 
	            dis(y2000,x2000,y103,x103)
                dis(y2000,x2000,y101,x101)
                dis(y3000,x3000,y100,x100)
                dis(y3000,x3000,y103,x103)
                dis(y100,x100,y1000,x1000)
                dis(y100,x100,y102,x102)
                dis(y100,x100,y2000,x2000)
                dis(y100,x100,y3000,x3000)
                dis(y101,x101,y103,x103)
                dis(y101,x101,y2000,x2000)
                dis(y102,x102,y1000,x1000)
                dis(y102,x102,y100,x100)
                dis(y102,x102,y2000,x2000)
                dis(y103,x103,y3000,x3000)
                dis(y103,x103,y2000,x2000)
                dis(y103,x103,y101,x101)
	            direction(y1000,x1000,y100,x100,w1000)
                direction(y1000,x1000,y102,x102,w1000) 
	            direction(y2000,x2000,y103,x103,w2000)
                direction(y2000,x2000,y101,x101,w2000)
                direction(y3000,x3000,y100,x100,w3000)
                direction(y3000,x3000,y103,x103,w3000)
                direction(y100,x100,y1000,x1000,w100)
                direction(y100,x100,y102,x102,w100)
                direction(y100,x100,y2000,x2000,w100)
                direction(y100,x100,y3000,x3000,w100)
                direction(y101,x101,y103,x103,w101)
                direction(y101,x101,y2000,x2000,w101)
                direction(y102,x102,y1000,x1000,w102)
                direction(y102,x102,y100,x100,w102)
                direction(y102,x102,y2000,x2000,w102)
                direction(y103,x103,y3000,x3000,w103)
                direction(y103,x103,y2000,x2000,w103)
                direction(y103,x103,y101,x101,w103)];               
	 
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
    while i>8
            X_0(i)=(X_0(i)*200/pi);
    end
   break
end    
X =["y100;" "x100;" "y101;" "x101;" "y102;" "x102;" "y103;" "x103;" "w1000;" "w2000;" "w3000;" "w100;" "w101;" "w102;" "w103;"]';
Results=[X,X_0];

max_len = max(length(L),length(L_hat));
for i = 1:max_len
    while i>18
        L_hat(i)=L_hat(i)*200/pi;
        L(i)=L(i)*200/pi;
        break
    end
end
adj_obs=[L,L_hat];

% Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

% VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

% Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));
for i=1:length(s_X)
    while i>8
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
for i=1:length(s_L_hat)
    while i>18
        s_L_hat(i)=s_L_hat(i)*200/pi;
        break
    end
end

std_adj_obs =[adj_obs,s_L_hat];
        
% Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

% VC matrix of residuals
S_vv = s_0^2*Q_vv;              

% Standard deviation of the residuals
s_v = sqrt(diag(S_vv));
for i=1:length(s_v)
    while i>18
        s_v(i)=s_v(i)*200/pi;
        v(i)=v(i)*200/pi;
        break
    end
end
std_resid =[v,s_v];



