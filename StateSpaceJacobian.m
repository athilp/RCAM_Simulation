function [F, G] = StateSpaceJacobian(Func, x0,u0)
N = length(x0);
F = zeros(N,N);

for i = 1:N
for j = 1:N
F(i,j) = ScalarGradient(Func, x0, u0, i, j, 0);
end
end

M = length(u0);
G = zeros(N, M);

for i = 1:N
for j = 1:M
G(i,j) = ScalarGradient(Func, x0, u0, i, j, 1);
end
end

end