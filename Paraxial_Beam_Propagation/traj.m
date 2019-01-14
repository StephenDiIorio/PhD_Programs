% calculates trajectory of single electron

% u0  is initial coordinates (x0, x0prime)
% d1,d2 distances source to lens and lens to screen

% u and ul are coords after and at lens
function [u,ul] = traj(u0, E0, K, r, R, d1, d2)

%first drift space
Dr1 = [1.0 d1
       0.0 1];

u1 = Dr1*u0;

% now focus
% u(1) is x
M = EfocMat(E0, u1(1), K, r, R);

ul = M*u1;

%final drift space
Dr2 = [1.0 d2
       0.0 1];

u = Dr2*ul;

end
