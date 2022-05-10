Z = 2;
cs = CrossSection;
ion_level = 1;
% nocc = cs.GetMaxOccShell(Z, ion_level);

opts = odeset('RelTol',2.22045e-14,'AbsTol',1e-32,'OutputFcn',@odeprog,'Events',@odeabort);

c = 299792458;
mec2 = 510998.95;
r_e = 2.8179403262E-15;
eps_0 = 8.8541878128E-12;
e_charge = 1.602176634E-19;
me = 9.1093837015E-31;
toJ = 1 / 6.24e18;


% n_0 = w_p^2 * me * eps_0 / e_charge^2;
n_0 = 2.16e23 * 1e6;
w_p = sqrt(n_0 * e_charge * e_charge / me / eps_0);
max_t = 6*pi;

u2 = 50;
% e_eng = 211662.69544841108;
gamma = sqrt(u2 + 1);
e_eng = u2 / (gamma + 1) * mec2;


%%
dens_mono_si = @(t,n) double(sqrt(2) .* ...
                             [( n(1)*n(2) * si_sigmavee_delta(Z,ion_level,  cs,e_eng)) + ( n(1)*n(3) * si_sigmavee_delta(Z,ion_level+1,cs,e_eng));...
                              (-n(1)*n(2) * si_sigmavee_delta(Z,ion_level,  cs,e_eng));...
                              ( n(1)*n(2) * si_sigmavee_delta(Z,ion_level,  cs,e_eng)) + (-n(1)*n(3) * si_sigmavee_delta(Z,ion_level+1,cs,e_eng));...
                              ( n(1)*n(3) * si_sigmavee_delta(Z,ion_level+1,cs,e_eng))]);
                          
dens_mono_si_v = @(t,n) double([( n(1)*n(2) * v_si_sigmavee_delta(Z,ion_level,  cs,u2)) + ( n(1)*n(3) * v_si_sigmavee_delta(Z,ion_level+1,cs,u2));...
                                (-n(1)*n(2) * v_si_sigmavee_delta(Z,ion_level,  cs,u2));...
                                ( n(1)*n(2) * v_si_sigmavee_delta(Z,ion_level,  cs,u2)) + (-n(1)*n(3) * v_si_sigmavee_delta(Z,ion_level+1,cs,u2));...
                                ( n(1)*n(3) * v_si_sigmavee_delta(Z,ion_level+1,cs,u2))]);

[t,y] = ode45(dens_mono_si_v,[0 max_t / w_p],[0.1 * n_0; n_0; 0; 0],opts);

figure(1)
subplot(2,2,1)

plot(t,y(:,1),'-o',t,y(:,2),'-o',t,y(:,3),'-o',t,y(:,4),'-o')
axis tight
set(gca,'FontSize',16)
xlabel('t (s)');
ylabel('\rho (m^{-3})');
title(strcat('SI Monoenergetic, E = ', num2str(e_eng), 'eV'))
legend('n_e','n_0','n_1','n_2','Location','best')


dens_mono = @(t,n) double((sqrt(2) / (4 * pi) * ...
                           c / (r_e * w_p)) .* ...
                          [( n(1)*n(2) * sigmavee_delta(Z,ion_level,  cs,e_eng/mec2,w_p)) + ( n(1)*n(3) * sigmavee_delta(Z,ion_level+1,cs,e_eng/mec2,w_p));...
                           (-n(1)*n(2) * sigmavee_delta(Z,ion_level,  cs,e_eng/mec2,w_p));...
                           ( n(1)*n(2) * sigmavee_delta(Z,ion_level,  cs,e_eng/mec2,w_p)) + (-n(1)*n(3) * sigmavee_delta(Z,ion_level+1,cs,e_eng/mec2,w_p));...
                           ( n(1)*n(3) * sigmavee_delta(Z,ion_level+1,cs,e_eng/mec2,w_p))]);

dens_mono_v = @(t,n) double((1 / (4 * pi) * ...
                             c / (r_e * w_p)) .* ...
                            [( n(1)*n(2) * v_sigmavee_delta(Z,ion_level,  cs,u2,w_p)) + ( n(1)*n(3) * v_sigmavee_delta(Z,ion_level+1,cs,u2,w_p));...
                             (-n(1)*n(2) * v_sigmavee_delta(Z,ion_level,  cs,u2,w_p));...
                             ( n(1)*n(2) * v_sigmavee_delta(Z,ion_level,  cs,u2,w_p)) + (-n(1)*n(3) * v_sigmavee_delta(Z,ion_level+1,cs,u2,w_p));...
                             ( n(1)*n(3) * v_sigmavee_delta(Z,ion_level+1,cs,u2,w_p))]);

[t,y] = ode45(dens_mono_v,[0 max_t],[0.1 * n_0 / n_0; n_0 / n_0; 0; 0],opts);

subplot(2,2,3)

plot(t,y(:,1),'-o',t,y(:,2),'-o',t,y(:,3),'-o',t,y(:,4),'-o')
axis tight
set(gca,'FontSize',16)
xlabel('\omega_p t');
ylabel('\rho / n_0');
title(strcat('Norm Monoenergetic, E = ', num2str(e_eng), 'eV'))
legend('n_e','n_0','n_1','n_2','Location','best')



si_val1 = double(si_sigmavee_max(Z,ion_level,cs,e_eng));
si_val2 = double(si_sigmavee_max(Z,ion_level+1,cs,e_eng));
dens_max_si = @(t,n) double((2 * sqrt(2/pi/me) * ...
                             (e_eng * toJ)^(-3/2) * ...
                             toJ^2) .* ...
                            [( n(1)*n(2) * si_val1) + ( n(1)*n(3) * si_val2);...
                             (-n(1)*n(2) * si_val1);...
                             ( n(1)*n(2) * si_val1) + (-n(1)*n(3) * si_val2);...
                             ( n(1)*n(3) * si_val2)]);

[t,y] = ode45(dens_max_si,[0 max_t / w_p],[0.1 * n_0; n_0; 0; 0],opts);

subplot(2,2,2)

plot(t,y(:,1),'-o',t,y(:,2),'-o',t,y(:,3),'-o',t,y(:,4),'-o')
axis tight
set(gca,'FontSize',16)
xlabel('t (s)');
ylabel('\rho (m^{-3})');
title(strcat('SI Maxwellian, k_BT = ', num2str(e_eng), 'eV'))
legend('n_e','n_0','n_1','n_2','Location','best')


val1 = double(sigmavee_max(Z,ion_level,cs,e_eng/mec2,w_p));
val2 = double(sigmavee_max(Z,ion_level+1,cs,e_eng/mec2,w_p));
dens_max = @(t,n) double((sqrt(2/pi)/2/pi *...
                          c / r_e / w_p) * ...
                         (e_eng/mec2)^(-3/2) .* ...
                         [( n(1)*n(2) * val1) + ( n(1)*n(3) * val2);...
                          (-n(1)*n(2) * val1);...
                          ( n(1)*n(2) * val1) + (-n(1)*n(3) * val2);...
                          ( n(1)*n(3) * val2)]);

[t,y] = ode45(dens_max,[0 max_t],[0.1 * n_0 / n_0; n_0 / n_0; 0; 0],opts);

subplot(2,2,4)

plot(t,y(:,1),'-o',t,y(:,2),'-o',t,y(:,3),'-o',t,y(:,4),'-o')
axis tight
set(gca,'FontSize',16)
xlabel('\omega_p t');
ylabel('\rho / n_0');
title(strcat('Norm Maxwellian, k_BT = ', num2str(e_eng), 'eV'))
legend('n_e','n_0','n_1','n_2','Location','best')


function val = sigmavee_delta(Z,ion_level,cs,e_eng,w_p)
    nocc = cs.GetMaxOccShell(Z, ion_level);
    
    val = 0;
    for ii = 1:nocc
        val = val + (cs.SingleCrossSecCalc_Norm(e_eng, Z, ion_level, ii, w_p) * sqrt(e_eng));
    end
end

function val = v_sigmavee_delta(Z,ion_level,cs,u2,w_p)
    nocc = cs.GetMaxOccShell(Z, ion_level);
    
    gamma = sqrt(u2 + 1);
    e_eng = u2 / (gamma + 1);
    
    val = 0;
    for ii = 1:nocc
        val = val + (cs.SingleCrossSecCalc_Norm(e_eng, Z, ion_level, ii, w_p) * sqrt(u2));
    end
end

function val = si_sigmavee_delta(Z,ion_level,cs,e_eng) % assumes e_eng comes in units of [eV]
    me = 9.1093837015E-31;
    toJ = 1 / 6.24e18;
    
    nocc = cs.GetMaxOccShell(Z, ion_level);
    
    val = 0;
    for ii = 1:nocc
        val = val + (cs.SingleCrossSecCalc(e_eng, Z, ion_level, ii) * sqrt(e_eng * toJ / me));
    end
end

function val = v_si_sigmavee_delta(Z,ion_level,cs,u2) % assumes u2 comes in normalized units
    mec2 = 510998.95;
    me = 9.1093837015E-31;
    c = 299792458;
    
    nocc = cs.GetMaxOccShell(Z, ion_level);
    
    u = sqrt(u2);   % convert to si units for calculation
    u = u * me * c; % this is momentum,
    v = u / me;
    
    gamma = sqrt(u2 + 1);
    e_eng = u2 / (gamma + 1) * mec2; % convert to si units [eV] for calculation
    
    val = 0;
    for ii = 1:nocc
        val = val + (cs.SingleCrossSecCalc(e_eng, Z, ion_level, ii) * v);
    end
end

function val = sigmavee_max(Z,ion_level,cs,kbT,w_p)
    syms T
    
    nocc = cs.GetMaxOccShell(Z, ion_level);
    
    val = 0;
    for ii = 1:nocc
        eqn = cs.GetEquation_Norm(Z,ion_level,ii,w_p);
        val = val + (vpaintegral(eqn * T * exp(-T/kbT),T,cs.GetBindingEnergy_Norm(Z, ion_level, ii),Inf, 'RelTol',1e-32, 'AbsTol',0));
    end
end

function val = si_sigmavee_max(Z,ion_level,cs,kbT)
    syms T
    
    nocc = cs.GetMaxOccShell(Z, ion_level);
    
    val = 0;
    for ii = 1:nocc
        eqn = cs.GetEquation(Z,ion_level,ii);
        val = val + (vpaintegral(eqn * T * exp(-T/kbT),T,cs.GetBindingEnergy(Z, ion_level, ii),Inf, 'RelTol',1e-32, 'AbsTol',0));
    end
end


%%

function [value,isterminal,direction]=odeabort(t,S,varargin)
%Other Events Set Here...ie:
% value(2)=max(abs(S(1:18)))-pi/2;
%Test to see if 'simstop' box is closed
    value(1) = double(ishandle(95));
    isterminal = 1;
    direction = 0;
end

function status = odeprog(t,y,flag,varargin)
%status = odebarplot(t,y,flag,varargin)
%   ODE progress display function with interrupt control
%   Displays a vertical bar plot that fills as the simulation
%   nears its completion.  Also displays time ellapsed and estimates
%   time remaining in the simulation.  To avoid computation burden
%   refreshes are limited to every 0.5 seconds.
%
% Tim Franklin
% Virginia Tech
% Jesse Norris
% Wake Forrest
% May 2006
    global odeprogglobvar
    if nargin < 3 || isempty(flag) 
        if(etime(clock,odeprogglobvar(8:13))>0.5)
            tfin = odeprogglobvar(1);
            sstrt = odeprogglobvar(2:7);
            figure(95); 
            perc = t(end) / tfin;
            area([t(end) tfin-t(end);t(end) tfin-t(end)]);
            title([num2str(perc*100) '%']);
            set(findobj('Tag','eltime'),'String',etimev(clock,sstrt));
            set(findobj('Tag','esttime'),'String',etimev(etime(clock,sstrt)/perc*(1-perc)));
            odeprogglobvar(8:13) = clock;
        end
    else
        switch(flag)
        case 'init'  
            odeprogglobvar = zeros(1,13);
            odeprogglobvar(1) = t(end);
            odeprogglobvar(2:7) = clock;
            odeprogglobvar(8:13) = clock;
            tfin = odeprogglobvar(1);
            sstrt = odeprogglobvar(2:7);
            figure(95); 
            set(gcf,'Position',[4,40,100,500]);
            axes('Position',[0.5,0.25,0.25,0.6]);
            axis([1,2,0,tfin]);
            set(gca,'XTickLabel',[],'NextPlot','replacechildren');
            ylabel('Simulation Progress - Time (s)');
            title('0%');
            area([0 tfin;0 tfin]);
            uicontrol('Style', 'pushbutton', 'String', 'Abort','Position', [7 460 90 30], 'Callback', 'close(gcf)')
            uicontrol('Style', 'text', 'String', 'Ellapsed Time','Position', [7 100 90 15])
            uicontrol('Style', 'text', 'Tag', 'eltime', 'String', etimev(clock,sstrt),'Position', [7 80 90 15])
            uicontrol('Style', 'text', 'String', 'Time Remaining','Position', [7 60 90 15])
            uicontrol('Style', 'text', 'Tag', 'esttime', 'String', num2str(inf),'Position', [7 40 90 15])
            pause(0.1);
        case 'done'    
            if(ishandle(95))
                close(95);
            end
        end
    end
    status = 0;
    drawnow;
end
function [S] = etimev(t1,t0)
%ETIMEV  Verbose Elapsed time.
%   ETIMEV(T1,T0) returns string of the days, hours, minutes, seconds that have elapsed 
%   between vectors T1 and T0.  The two vectors must be six elements long, in
%   the format returned by CLOCK:
%
%       T = [Year Month Day Hour Minute Second]
%   OR
%   ETIMEV(t), t in seconds
    if(exist('t1') && exist('t0') && length(t1)>2 && length(t0)>2)
        t=etime(t1,t0);     %Time in seconds
        if(t<0)
            t = -t;
        end
    elseif(length(t1)==1)
        t = t1;
    else
        t=0;
    end
    days = floor(t/(24*60*60));
    t = t - days*24*60*60;
    hours = floor(t/(60*60));
    t = t - hours*60*60;
    mins = floor(t/60);
    t = floor(t-mins*60);
    if(days>0)
        S=[num2str(days) 'd ' num2str(hours) 'h ' num2str(mins) 'm ' num2str(t) 's'];
    elseif(hours>0)
        S=[num2str(hours) 'h ' num2str(mins) 'm ' num2str(t) 's'];
    elseif(mins>0)
        S=[num2str(mins) 'm ' num2str(t) 's'];
    else
        S=[num2str(t) 's'];
    end
end