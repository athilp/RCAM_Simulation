function F0 = cost_trim(Z)
%-----------------------------------------------------------------------------------------------
% Brief: Computes the cost function as a function of state and input
% motion.
% Details:
%    None
% 
% Syntax:  
%     f0 = cost_trim(z)
% 
% Inputs:
%    Z - [15,1] size, [double] type, State and Input parameters
% 
% Outputs:
%    f0 - [1,1] size,[double] type, Cost 
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
Vd = 85;

X = Z(1:9);
U = Z(10:14);

if U(4) < 0
    F0 = 10^6;
    return;
end
if U(5) < 0
    F0 = 10^6;
    return
end

xdot = RCAM_Model(X,U);
theta = X(8);
Va = sqrt(X(1)^2+X(2)^2+X(3)^2);
alpha = atan2(X(3), X(1));
gammma = theta - alpha;

C1v = [cos(X(9)) sin(X(9)) 0; -sin(X(9)) cos(X(9)) 0; 0 0 1];
C21 = [cos(X(8)) 0 -sin(X(8)); 0 1 0; sin(X(8)) 0 cos(X(8))];
Cb2 = [1 0 0; 0 cos(X(7)) sin(X(7)); 0 -sin(X(7)) cos(X(7))];
Cbv = Cb2 * C21 * C1v;
V_v = Cbv' * [X(1);X(2);X(3)];
hdot = -V_v(3);

Q = [xdot(4:6);hdot;Va-Vd;gammma;X(7)-0.26179939;X(9)];

%Exclude x10 state in the cost function
H = diag(ones(1, 8));
F0 = Q'*H*Q;
end