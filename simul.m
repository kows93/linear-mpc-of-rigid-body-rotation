clc; clear; close all;
h = 0.1;
N = 4;
Simul_time = 10;
figure_option = 1; % {1: Desired & Current trajectories, 2: tracking error}

x0 = [-0.5;0.3;0.2;0.05;0;1.15];
x_history = x0';
simul_t = 0;

for kh = 0:h:Simul_time-h
    xd = reference(kh);
    z = x0-xd(1:6);
    v = mpc(kh,z,h,N);
    [t,x]=ode45(@(t,x) sys(t,x,v,kh),[kh kh+h], x0);
    x0 = x(end,1:6)';
    x_history = [x_history; x(2:end,1:6)];
    simul_t = [simul_t;t(2:end)];
end
    

desired = zeros(size(x_history,1),6);

for i = 1:length(simul_t)
    [x_d,~] = reference(simul_t(i));
    desired(i,:) = x_d(1:6)';
end

if figure_option == 1
    figure(1);
    plot(simul_t,desired(:,1),simul_t,x_history(:,1)); grid on;
    legend('\phi_d','\phi');
    xlim([0 10]);

    figure(2);
    plot(simul_t,desired(:,2),simul_t,x_history(:,2)); grid on;
    legend('\theta_d','\theta');
    xlim([0 10]);

    figure(3);
    plot(simul_t,desired(:,3),simul_t,x_history(:,3)); grid on;
    legend('\psi_d','\psi');
    xlim([0 10]);

    figure(4);
    plot(simul_t,desired(:,4),simul_t,x_history(:,4)); grid on;
    legend('\Omega_{1,d}','\Omega_1');
    xlim([0 10]);

    figure(5);
    plot(simul_t,desired(:,5),simul_t,x_history(:,5)); grid on;
    legend('\Omega_{2,d}','\Omega_2');
    xlim([0 10]);

    figure(6);
    plot(simul_t,desired(:,6),simul_t,x_history(:,6)); grid on;
    legend('\Omega_{3,d}','\Omega_3');
    xlim([0 10]);

elseif figure_option == 2
    figure(1);
    plot(simul_t,abs(x_history(:,1)-desired(:,1))); grid on;
    ylabel('|\phi - \phi_d|');
    xlim([0 10]);

    figure(2);
    plot(simul_t,abs(x_history(:,2)-desired(:,2))); grid on;
    ylabel('|\theta - \theta_d|');
    xlim([0 10]);

    figure(3);
    plot(simul_t,abs(x_history(:,3)-desired(:,3))); grid on;
    ylabel('|\psi - \psi_d|');
    xlim([0 10]);

    figure(4);
    plot(simul_t,abs(x_history(:,4)-desired(:,4))); grid on;
    ylabel('|\Omega_1 - \Omega_{1,d}|');
    xlim([0 10]);

    figure(5);
    plot(simul_t,abs(x_history(:,5)-desired(:,5))); grid on;
    ylabel('|\Omega_2 - \Omega_{2,d}|');
    xlim([0 10]);

    figure(6);
    plot(simul_t,abs(x_history(:,6)-desired(:,6))); grid on;
    ylabel('|\Omega_3 - \Omega_{3,d}|');
    xlim([0 10]);
end