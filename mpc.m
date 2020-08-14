function v = mpc(t,z,h,N)
Q = diag([100 100 100 10 10 10]);
Q = kron(eye(N-1),Q);
R = 0.1*eye(3*N);
S = diag([100 100 100 10 10 10]);
[A_eq, b_eq] = dynamics_const(t,z,h,N);
z0 = initialize(t,z,h,N);
v0 = zeros(3*N,1);
cost = @(x)x(1:6*(N-1))'*Q*x(1:6*(N-1)) + x(6*(N-1)+1:6*N)'*S*x(6*(N-1)+1:6*N) + x(6*N+1:end)'*R*x(6*N+1:end);
[x, ~]=fmincon(cost, [z0;v0], [],[],A_eq, b_eq);
v = x(6*N+1:6*N+1+2);
end
