function [x_d,u_d] = reference(t)
I = diag([4.250 4.337 3.664]);
e_d = [t/10; t/10; t/10];
e_dt = [1/10; 1/10; 1/10];
e_ddt = [0; 0; 0];

S = [1 0 -sin(e_d(2)); 0 cos(e_d(1)) cos(e_d(2))*sin(e_d(1)); 0 -sin(e_d(1)) cos(e_d(2))*cos(e_d(1))];
o_d = S*e_dt;

Sd = [0 0 -cos(e_d(2))*e_dt(2); 0 -sin(e_d(1))*e_dt(1) -sin(e_d(2))*sin(e_d(1))*e_dt(2) + cos(e_d(2))*cos(e_d(1))*e_dt(1);
    0 -cos(e_d(1))*e_dt(1) -sin(e_d(2))*cos(e_d(1))*e_dt(2) - cos(e_d(2))*sin(e_d(1))*e_dt(1)];
o_dt = Sd*e_dt + S*e_ddt;

x_d = [e_d;o_d;e_dt;o_dt];
u_d = I*o_dt - cross((I*o_d),o_d);
end

