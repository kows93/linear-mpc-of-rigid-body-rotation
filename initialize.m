function z0 = initialize(t,z,T,N)
I = [4.250 4.337 3.664];
z_temp = z;
z0 = [];
for i=1:N
    [x_d,u_d] = reference(t+(i-1)*T);
    e_d = x_d(1:3);
    o_d = x_d(4:6);    
    A_temp = [cos(e_d(1))*tan(e_d(2))*o_d(2)-sin(e_d(1))*tan(e_d(2))*o_d(3) sin(e_d(1))*sec(e_d(2))^2*o_d(2)+cos(e_d(1))*sec(e_d(2))^2*o_d(3) 0 1 sin(e_d(1))*tan(e_d(2)) cos(e_d(1))*tan(e_d(2));
        -sin(e_d(1))*o_d(2)-cos(e_d(1))*o_d(3) 0 0 0 cos(e_d(1)) -sin(e_d(1));
        cos(e_d(1))/cos(e_d(2))*o_d(2)-sin(e_d(1))/cos(e_d(2))*o_d(3) sin(e_d(1))*sec(e_d(2))*tan(e_d(2))*o_d(2)+cos(e_d(1))*sec(e_d(2))*tan(e_d(2))*o_d(3) 0 0 sin(e_d(1))/cos(e_d(2)) cos(e_d(1))/cos(e_d(2));
        0 0 0 0 (I(2)-I(3))/I(1)*o_d(3) (I(2)-I(3))/I(1)*o_d(2);
        0 0 0 (I(3)-I(1))/I(2)*o_d(3) 0 (I(3)-I(1))/I(2)*o_d(1);
        0 0 0 (I(1)-I(2))/I(3)*o_d(2) (I(1)-I(2))/I(3)*o_d(1) 0];
        
    W_temp = ([[1 sin(e_d(1))*tan(e_d(2)) cos(e_d(1))*tan(e_d(2)); 0 cos(e_d(1)) -sin(e_d(1));
0 sin(e_d(1))/cos(e_d(2)) cos(e_d(1))/cos(e_d(2))]*o_d;
    pinv(diag(I))*(cross((diag(I)*o_d),o_d) + u_d)]-x_d(7:end))*T;
    A_temp = A_temp*T+eye(6);
    z_temp = A_temp*z_temp + W_temp;
    z0 = [z0;z_temp];
end
end
    