function [XDOT] = RCAM_Model(X, U)
%-----------------------------------------------------------------------------------------------
% Brief: Computes the derivatives of the 6 DOF aircraft equations of
% motion.
% Details:
%    None
% 
% Syntax:  
%     xdot(1:10) = RCAM_Model(X, U)
% 
% Inputs:
%    X - [10,1] size, [double] type,State Vector
%    U - [5,1] size,[double] type,Input Vector
% 
% Outputs:
%    xdot - [10,1] size,[double] type, Derivatives for ODE45 or Simulink
%    Integration
% 
% Example: 
%    None
% 
% See also: None

% Author:                          Athil George
% Email:                           athilg@utexas.edu
% Created:                         12-Oct-2023
% Version history revision notes:
%                                  None
%-----------------------------------------------------------------------------------------------
%% 1. State and Control Defintions

%Define the states
x1 = X(1);
x2 = X(2);
x3 = X(3);
x4 = X(4);
x5 = X(5);
x6 = X(6);
x7 = X(7);
x8 = X(8);
x9 = X(9);
% x10 = X(10);

%Define the controls
u1 = U(1);
u2 = U(2);
u3 = U(3);
u4 = U(4);
u5 = U(5);

%-------------------------------------------------------------------------------------------
%% 2. Key Parameters

%Define vehicle parameters
m = 120000;
cbar = 6.6;
lt = 24.8;
S = 260;
St = 64;

%Center of gravity locations
Xcg = 0.23*cbar;
Ycg = 0;
Zcg = 0.10*cbar;

%Aerodynamic Center Locations
Xac = 0.12*cbar;
Yac = 0;
Zac = 0;

%Engine Applciation locations
Xapt1 = 0;
Yapt1 = -7.94;
Zapt1 = -1.9;

Xapt2 = 0;
Yapt2 = 7.94;
Zapt2 = -1.9;

% ---------------------------------------------------------------------------------------
%% 3. Air density Model

% if x10 > 1000000
%     rho = 3.019e-15;
% elseif x10 > 900000
%     rho = 5.245e-15;
% elseif x10 > 800000
%     rho = 1.170e-14;
% elseif x10 > 700000
%     rho = 3.614e-14;
% elseif x10 > 600000
%     rho = 1.454e-13;
% elseif x10 > 500000
%     rho = 6.967e-13;
% elseif x10 > 450000
%     rho = 1.585e-12;
% elseif x10 > 400000
%     rho = 3.725e-12;
% elseif x10 > 350000
%     rho = 9.518e-12;
% elseif x10 > 300000
%     rho = 2.418e-11;
% elseif x10 > 250000
%     rho = 7.248e-11;
% elseif x10 > 200000
%     rho = 2.789e-10;
% elseif x10 > 180000
%     rho = 5.464e-10;
% elseif x10 > 150000
%     rho = 2.070e-9;
% elseif x10 > 140000
%     rho = 3.845e-9;
% elseif x10 > 130000
%     rho = 8.484e-9;
% elseif x10 > 120000
%     rho = 2.438e-8;
% elseif x10 > 110000
%     rho = 9.661e-8;
% elseif x10 > 100000
%     rho = 5.297e-7;
% elseif x10 > 90000
%     rho = 3.396e-6;
% elseif x10 > 80000
%     rho = 1.905e-5;
% elseif x10 > 70000
%     rho = 8.770e-5;
% elseif x10 > 60000
%     rho = 3.206e-4;
% elseif x10 > 50000
%     rho = 1.057e-3;
% elseif x10 > 40000
%     rho = 3.972e-3;
% elseif x10 > 30000
%     rho = 1.774e-2;
% elseif x10 > 25000
%     rho = 3.899e-2;
% else
%     rho = 1.225;
% end
rho = 1.225;
% ---------------------------------------------------------------------------------------

%Define aerodynamic parameters
g = 9.81;
depsda = 0.25;
alpha_L0 = -11.5 * pi/180;
n = 5.5;
a3 = -768.5;
a2 = 609.2;
a1 = -155.2;
a0 = 15.212;
alpha_switch = 14.5 * pi/180;

%Define the airspeed
Va = sqrt(x1^2 + x2^2 + x3^2);

%Get the stability angles
alpha = atan2(x3,x1);
beta = asin(x2/Va);

%Dynamic Pressure
Q = 0.5*rho*Va^2;

%Define the angular velocity and velocity vector in the body frame
wbe_b = [x4;x5;x6];
V_b = [x1;x2;x3];

%Define the lift curve slope
if alpha <= alpha_switch
    CL_wb = n*(alpha - alpha_L0);
else
    CL_wb = a3*alpha^3 + a2*alpha^2 + a1*alpha + a0;
end

%Define the downwash angles due to angle of attakc and lift for the tail
epsilon = depsda*(alpha - alpha_L0);
alpha_t = alpha - epsilon + u2 + 1.3*x5*lt/Va;
CL_t = 3.1*(St/S)*alpha_t;

%Get the forces in the stability axes
CL = CL_wb + CL_t;
CD = 0.13+ 0.07*(5.5*alpha+0.654)^2;
CY = -1.6*beta + 0.24*u3;
FA_s = [-CD*Q*S; CY*Q*S;-CL*Q*S];

%Transform into the body axes
C_bs = [cos(alpha) 0 -sin(alpha); 0 1 0; sin(alpha) 0 cos(alpha)];
FA_b = C_bs*FA_s;

%Define the intercept moment contribution
eta11 = -1.4*beta;
eta21 = -0.59 - (3.1*(St*lt)/(S*cbar))*(alpha - epsilon);
eta31 = (1 - alpha *(180/(15*pi)))*beta;
eta = [eta11;eta21;eta31];

%Define the control and state moment authorities
dCMdx = (cbar/Va) * [-11 0 5; 0 (-4.03*(St*lt^2)/(S*cbar^2)) 0; 1.7 0 -11.5];
dCMdu = [-0.6 0 0.22; 0 (-3.1*(St*lt)/(S*cbar)) 0; 0 0 -0.63];

%Compute the moment coefficients in the body frame
CMac_b = eta + dCMdx * wbe_b + dCMdu * [u1;u2;u3];

%Compute the moment about the AC
MAac_b = CMac_b * Q * S * cbar;

%Define the vectors from an arbitrary origin to the CG and AC
rcg_b = [Xcg;Ycg;Zcg];
rac_b = [Xac;Yac;Zac];

%Moment transfer to the CG
MAcg_b = MAac_b + cross(FA_b, rcg_b - rac_b);

% ---------------------------------------------------------------------------------------
%% 4. Propulsion Effects

%Compute the forces from each engine given the throttle setting
F1 = u4*m*g;
F2 = u5*m*g;

a = 1.1772e-08;
N1 = (F1/a)^(-3.5);
N2 = (F1/a)^(-3.5);

N1_radps = 0.10472*N1;
N2_radps = 0.10472*N2;

br = 0.15+0.15+0.3+0.1;
rE = 1;
mE = 19504;
IE = 0.5 * (br*mE) * rE^2;

hrot = [IE * N1_radps; 0; 0];

FE1_b = [F1;0;0];
FE2_b = [F2;0;0];
FE_b = FE1_b + FE2_b;

%Position vector from ENGAPT to CG
mew1 = [Xcg - Xapt1; Yapt1 - Ycg; Zcg - Zapt1];
mew2 = [Xcg - Xapt2; Yapt2 - Ycg; Zcg - Zapt2];

%COmpute the moment about the center of gravity due to engine effects
MEcg1_b = cross(mew1, FE1_b);
MEcg2_b = cross(mew2, FE2_b);
MEcg_b = MEcg1_b + MEcg2_b;

% ---------------------------------------------------------------------------------------
%% 5. Gravity Effects

%Define the gravity vector in the body frame. This is rotated by a DCM
g_b = [-g*sin(x8); g*cos(x8)*sin(x7); g*cos(x8)*cos(x7)];

%Compute the Force
Fg_b = m*g_b;

%COmpute the i%nertial matri%x and its inverse
Ib = m*[50.07 0 -2.0923; 0 64 0; -2.0923 0 99.92];
invIb = (1/m)*[0.0249836 0 0.000523151; 0 0.015625 0; 0.000523151 0 0.010019];

%Total foroce due to aerodynamic, gravity, and propulsion
F_b = Fg_b + FE_b + FA_b;

x1to3dot = (1/m)*F_b - cross(wbe_b, V_b);

%Total moment due to aerodynamic, gravity, and propulsion
Mcg_b = MAcg_b + MEcg_b;

wbe_b_CPE = [0 -wbe_b(3) wbe_b(2); wbe_b(3) 0 -wbe_b(1); -wbe_b(2) wbe_b(1) 0];
% x4to6dot = invIb*(Mcg_b - cross(wbe_b,Ib*wbe_b));
x4to6dot = invIb*(Mcg_b - wbe_b_CPE*(Ib*wbe_b + hrot));

%Use Poissions Kinematic Equations to get phi, theta, psi
H_phi = [1 sin(x7)*tan(x8) cos(x7)*tan(x8);0 cos(x7) -sin(x7);0 sin(x7)/cos(x8) cos(x7)/cos(x8)];
x7to9dot = H_phi*wbe_b;
% 
% %DCM -> NED frame to get rate of climb/descent
% C1v = [cos(x9) sin(x9) 0; -sin(x9) cos(x9) 0; 0 0 1];
% C21 = [cos(x8) 0 -sin(x8); 0 1 0; sin(x8) 0 cos(x8)];
% Cb2 = [1 0 0; 0 cos(x7) sin(x7); 0 -sin(x7) cos(x7)];
% Cbv = Cb2 * C21 * C1v;
% V_v = Cbv' * V_b;
% xdot10 = -V_v(3);

XDOT = [x1to3dot;x4to6dot;x7to9dot];
end






