% this makes time dependent map of E probe
% using thin lens approximation
function Efieldprobe

% distance particle source to lens
d1 = 2000;

% distance  lens to screen
d2 = 180*8000;

% number of particles
N = 10000;

%width of distribution in x
dx = 500.0;

% width in xprime (i.e. angle spread in rad)
dxprime = 0.0002;

K0 = 90/511; % peak kinetic energy of particles

delta = 10.0; % thickness of E field shell

E0 = 0.00015; 0.002; % Efield magnitude at
r0 = 50.0; % radius of shell

v_expand = 0.1; % "expansion velocity"

t = -50.0; % initial time
tmax = 200; % final time
dt = 1.0; % step size

showtraj = false; % show 50 trajectories
Etraj = 0.0000015;
rtraj = 100;

x = linspace(-500,500,100); % coordinate on screen

screen_t = [];
time_t = [];
E_t = [];
K_t = [];

if showtraj
    makeprofile(Etraj,K0,rtraj,delta,x,showtraj,dxprime,dx,N,d1,d2);
else
    while t<tmax
        % use a field of the form E = E0*r0/r
        % and r proportional to time

        r = v_expand*sqrt(t)*1000;

        if t>0
            E = E0*r0./r;
        else
            E=0.0;
        end

        K = (1-(t/70-1)*2.0/9.0)*K0;
        profile = makeprofile(E,K,r,delta,x,showtraj,dxprime,dx,N,d1,d2);

        E_t = [E_t (E.*(abs(x)>r).*(abs(x)<(r+delta)))'];
        screen_t = [screen_t profile'];
        K_t = [K_t K];
        time_t = [time_t t];
        t = t+dt
    end

    % subplot(211);
    subplot(111);
    imagesc(time_t,x*0.8/2/pi,screen_t);
    c = colorbar;
    c.FontSize = 16;
    xlabel('t (ps)', 'Fontsize', 18);
    ylabel('x ({\mu}m)', 'Fontsize', 18);
    title('Electron Density on Screen', 'Fontsize', 20);
    % subplot(212);
    % imagesc(time_t,x*0.8/2/pi,E_t);
    % c=colorbar;
    % c.FontSize = 16;
    % xlabel('time [ps]', 'Fontsize', 15);
    % ylabel('x [{\mu}m]', 'Fontsize', 15);
    % title('Electric field model', 'Fontsize', 20);
    figure;
    profile = sum(screen_t);
    subplot(111)
    plot(time_t,profile)
    % hAx = plotyy(time_t,profile,time_t,K_t);
    xlabel('t (ps)', 'Fontsize', 18);
    ylabel('Electron density on screen', 'Fontsize', 18);
    title('Lineout Along Streaking Axis', 'Fontsize', 20);
    % ylabel(hAx(1),'Electron density on screen', 'Fontsize', 15);
    % ylabel(hAx(2),'Electron energy', 'Fontsize', 15);

end

end
