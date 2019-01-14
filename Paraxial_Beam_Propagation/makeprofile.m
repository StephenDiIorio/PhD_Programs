% this generates a random distribution of particles
% calculates their trajectories through a thin shell field
% of radius R and thickness delta
% field strength E0 and kinetic energy K
% screen_x is the x_position grid on the screen
% ifplottraj is bool - whether to plot all the trajectories
% NOT VECTORIZED, SLOW - I know!

function profile = makeprofile(E0, K, R, delta, screen_x, ifplottraj, dxprime, dx, Npar, d1, d2)

x0list = [];
N = Npar;
if ifplottraj
    N = 50;
end

for ii = 1:N
    x0 = dx*randn;
    x0prime = dxprime*randn;
    u0 = [x0;x0prime];
    [u,ul] = traj(u0, E0, K, R, R+delta, d1, d2);

    % display trajectories
    if (ifplottraj)
        plot([-d1 0.0 d2],[u0(1) ul(1) u(1)]);
        ylim([-dx dx]);
        hold on
    end

   % add to list of x0s
   x0list = [x0list u(1)];
end

dxs = screen_x(2) - screen_x(1);
[N_w,jj_l,qj_l] = weighting(N, x0list - min(screen_x), dxs);

Nx = length(screen_x);
profile = zeros(1,Nx);
for ll=1:N_w
    for ii=1:length(x0list)
        if (jj_l(ll,ii) >0 && jj_l(ll,ii) <=Nx)
            profile(jj_l(ll,ii)) = profile(jj_l(ll,ii)) + qj_l(ll,ii);
        end
    end

    %plot(screen_x,profile);
end

end
