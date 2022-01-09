function [dot_x] = my_rightSide(x, params)
dot_x = zeros(13, 1);
R_vect = x(1:3);
V_vect = x(4:6);
w_vect_IF = x(7:9);
Q = x(10:13)';
IF2BF = quat2dcm(Q)';

j3_IF = R_vect / norm(R_vect);
j3_BF = IF2BF * j3_IF;
mu = params.mu;
M_vect_BF = 3*mu/norm(R_vect).^3 * cross(j3_BF, params.J * j3_BF);
w_vect_BF = IF2BF * w_vect_IF;
dot_w_vect = params.inv_J*(M_vect_BF - cross(w_vect_BF, params.J * w_vect_BF));

dot_x(1:3) = V_vect;
dot_x(4:6) = -mu/norm(R_vect).^3 * R_vect;
dot_x(7:9) = IF2BF' * dot_w_vect;
dot_x(10:13) = 0.5*quatmultiply(Q, [0, w_vect_IF']);
end
