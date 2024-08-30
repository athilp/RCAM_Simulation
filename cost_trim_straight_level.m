function F0 = cost_trim_straight_level(Z)
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
Vd = 230;

X = Z(1:9);
U = Z(10:14);

xdot = RCAM_Model(X,U);
theta = X(8);
Va = sqrt(X(1)^2+X(2)^2+X(3)^2);
alpha = atan2(X(3), X(1));
gammma = theta - alpha;

Q = [xdot(1:9);Va-Vd;gammma;X(2);X(7);X(9)];

%Exclude x10 state in the cost function
H = diag(ones(1, 14));
F0 = Q'*H*Q;
end