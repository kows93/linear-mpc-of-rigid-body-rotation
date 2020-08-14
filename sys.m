function dx = sys(t,x,v,kh)
e = x(1:3); 
o = x(4:6); 

[~,u_d] = reference(kh); 
tau = v+u_d;

Q = [1 sin(e(1))*tan(e(2)) cos(e(1))*tan(e(2)); 0 cos(e(1)) -sin(e(1)); 0 sin(e(1))/cos(e(2)) cos(e(1))/cos(e(2))];
I = diag([4.250 4.337 3.664]); 
de = Q*o;
do = pinv(I)*(cross((I*o),o) + tau);

dx = [de; do];
end