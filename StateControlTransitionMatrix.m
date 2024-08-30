function [PHI, GAMMA] = StateControlTransitionMatrix(F, G, dt)
n = size(F,1);
PHI = expm(F * dt);
GAMMA = (eye(n,n) - PHI) * G;
end