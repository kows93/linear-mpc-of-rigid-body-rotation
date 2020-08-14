function [A_eq, b_eq] = dynamics_const(t,z,h,N)
A = zeros(6*N);
B = zeros(6*N,3*N);
W = [];
I = [4.250 4.337 3.664];
for i=0:N-1
    [x_d,u_d] = reference(t+(i)*h);
    e_d = x_d(1:3);
    o_d = x_d(4:6);    
    A_temp = [cos(e_d(1))*tan(e_d(2))*o_d(2)-sin(e_d(1))*tan(e_d(2))*o_d(3) sin(e_d(1))*sec(e_d(2))^2*o_d(2)+cos(e_d(1))*sec(e_d(2))^2*o_d(3) 0 1 sin(e_d(1))*tan(e_d(2)) cos(e_d(1))*tan(e_d(2));
        -sin(e_d(1))*o_d(2)-cos(e_d(1))*o_d(3) 0 0 0 cos(e_d(1)) -sin(e_d(1));
        cos(e_d(1))/cos(e_d(2))*o_d(2)-sin(e_d(1))/cos(e_d(2))*o_d(3) sin(e_d(1))*sec(e_d(2))*tan(e_d(2))*o_d(2)+cos(e_d(1))*sec(e_d(2))*tan(e_d(2))*o_d(3) 0 0 sin(e_d(1))/cos(e_d(2)) cos(e_d(1))/cos(e_d(2));
        0 0 0 0 (I(2)-I(3))/I(1)*o_d(3) (I(2)-I(3))/I(1)*o_d(2);
        0 0 0 (I(3)-I(1))/I(2)*o_d(3) 0 (I(3)-I(1))/I(2)*o_d(1);
        0 0 0 (I(1)-I(2))/I(3)*o_d(2) (I(1)-I(2))/I(3)*o_d(1) 0];
    A_temp = A_temp*h+eye(6);
    if i == 0
        Ak = A_temp;
    else
        A(6*(i)+1:6*(i+1),6*(i-1)+1:6*(i)) = A_temp;        
    end
    
    A(6*(i)+1:6*(i+1), 6*(i)+1:6*(i+1)) = -eye(6);
    
    B_temp = [zeros(3,3);pinv(diag(I))];
    B(6*(i)+1:6*(i+1),3*(i)+1:3*(i+1)) = B_temp*h;
    
    W_temp=([[1 sin(e_d(1))*tan(e_d(2)) cos(e_d(1))*tan(e_d(2)); 0 cos(e_d(1)) -sin(e_d(1));
        0 sin(e_d(1))/cos(e_d(2)) cos(e_d(1))/cos(e_d(2))]*o_d;
        pinv(diag(I))*(cross((diag(I)*o_d),o_d) + u_d)]-x_d(7:end))*h;
    if i == 0
        W_temp = W_temp + Ak*z;
    end
    W = [W; W_temp];
end

A_eq = [A B];
b_eq = -W;
end


