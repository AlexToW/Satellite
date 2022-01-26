clear all;
close all;
clc;

dt = 0.5;
t_end = 15000;
t = 0:dt:t_end;
N = length(t);
params = struct();
A = 2;
B = 2;
C = 4;
params.J = diag([A,B,C]);
params.inv_J = diag([1/A, 1/B, 1/C]);
params.m = 1000;
params.mu = 3.987*10.^(14);
dot_theta0 = 0.1;
dot_psi0 = 0.3;
dot_phi0 = 0.1;

R_0 = 6.371 * 10.^6; % радиус Земли
alpha_0 = pi/6;
R_0_vect = R_0 * [0; cos(alpha_0); sin(alpha_0)];
V_0 = sqrt(params.mu / R_0); % первая космическая
V_0_vect = V_0 * [0;-sin(alpha_0);cos(alpha_0)];

psi_0 = 0;
theta_0 = 0;
phi_0 = 0;
B_psi_0 = [cos(psi_0), -sin(psi_0), 0;
         sin(psi_0), cos(psi_0), 0;
         0,         0,       1];
B_theta_0 = [1,       0,        0;
            0,    cos(theta_0), -sin(theta_0);
            0,    sin(theta_0), cos(theta_0)];
B_phi_0 = [cos(phi_0), -sin(phi_0), 0;
         sin(phi_0), cos(phi_0),  0;
         0,         0,        1];
B_0 = B_psi_0 * B_theta_0 * B_phi_0;
% орбитальная СК
j1_0 = V_0_vect/norm(V_0_vect);
j3_0 = R_0_vect/norm(R_0_vect);
j2_0 = cross(j3_0, j1_0);
j2_0 = j2_0/norm(j2_0);
j1_0 = cross(j2_0, j3_0);
A_0 = [j1_0, j2_0, j3_0];
C_0 = A_0 * B_0;


% w_0_vect -- абсолютная угловая скорость в ССК при t=0
w_0_vect_rel = [1;0;0]; % относительная скорость в ССК при t=0
w_0_vect = w_0_vect_rel + C_0 * cross(R_0_vect, V_0_vect)/norm(R_0_vect).^2;
Q_0 = dcm2quat(C_0');

% фазовый вектор
x = zeros(13, N);
x(:, 1) = [R_0_vect; V_0_vect; w_0_vect; Q_0'];

for i = 1:N-1
    k1 = my_rightSide(x(:, i), params);
    k2 = my_rightSide(x(:, i) + k1*dt/2, params);
    k3 = my_rightSide(x(:, i) + k2*dt/2, params);
    k4 = my_rightSide(x(:, i) + k3*dt, params);
    x(:, i+1) = x(:, i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    x(10:13, i+1) = x(10:13, i+1) / norm(x(10:13, i+1));
end


Omegas_orb = zeros(3, N); % угловая скорость орбиталки в ИСК
omegas_relative = zeros(3, N); % относительная угл. ск-ть в ССК
Omegas_orb(1:3, 1) = cross(R_0_vect, V_0_vect)/norm(R_0_vect).^2;
for i=2:N
    V_vect = x(4:6, i); % в ИСК
    R_vect = x(1:3, i); % в ИСК
    C_tmp = quat2dcm(x(10:13, i)')';
    Omegas_orb(1:3, i) = cross(R_vect, V_vect)/norm(R_vect).^2;
    w_abs_BF = x(7:9, i);
    omegas_relative(1:3, i) = w_abs_BF - C_tmp * Omegas_orb(1:3, i);
    fprintf("(");
    for j=1:3
        fprintf("%f, ", omegas_relative(j, i));
    end
    fprintf(")\n");
end

figure
hold on
grid on
axis equal
plot(x(2,:), x(3,:));




