% Weights particles to/from grid
function [N_w,jj_l,qj_l] = weighting(Nx, xi, dx)
% Second order weighting!
% particle shape is triangle of base 2dx
N_w = 3;

%First generate vector rounded xi/dx;
jj = floor(xi/dx+0.5);

Delta = (xi-jj*dx)/dx;

% find weight at j+1
qjp = 0.5*(0.5+Delta).^2;
% find weight at j-1
qjm = 0.5*(0.5-Delta).^2;

% and at j
qj = 1.0-(qjp+qjm);% (1.0*(jj+1)*dx-xi);

jjp = jj+1;
jjm = jj-1;

% apply boundary condition
% jj=jj+((jj<0) - (jj>Nx-1))*Nx;
% jjp=jjp+((jjp<0) - (jjp>Nx-1))*Nx;
% jjm=jjm+((jjm<0) - (jjm>Nx-1))*Nx;

% because of matlab indexing
jj = jj+1;
jjp = jjp+1;
jjm = jjm+1;

jj_l(1,:) = jjm;
jj_l(2,:) = jj;
jj_l(3,:) = jjp;

qj_l(1,:) = qjm;
qj_l(2,:) = qj;
qj_l(3,:) = qjp;

end
