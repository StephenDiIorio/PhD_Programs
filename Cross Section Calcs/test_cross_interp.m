%% Cross Sectioz Set Up

c = 299792458;
mec2 = 510998.95;
r_e = 2.8179403262E-15;
eps_0 = 8.8541878128E-12;
e_charge = 1.602176634E-19;
me = 9.1093837015E-31;
% n_0 = 2.16e23 * 1e6;
n_0 = 1.75e21 * 1e6;
w_p = sqrt(n_0 * e_charge * e_charge / me / eps_0);
toJ = 1 / 6.24e18;

Z = 13;
ion_level = 1;

cs = CrossSection;
num_points = 250;
e_max = 1e12 / mec2;
e_min = (cs.GetMinBindingEnergy(Z, ion_level) + 0.001) / mec2;
lin_slope = (num_points - 1) / (e_max - e_min);
log_slope = (num_points - 1) / log10(e_max / e_min);

occ_idx = cs.GetMaxOccShell(Z, ion_level);


%% Cross Section Test
% lin_energies = linspace(e_min, e_max, num_points);
energies = logspace(log10(e_min),log10(e_max),num_points);
sig = zeros(size(energies));
lost = zeros(size(energies));
transf = zeros(size(energies));
shell_sig = zeros(length(energies),occ_idx);


% test_e = 0.41421;
test_e = 1.17406480E-05;
num = linspace(1,num_points,num_points);
y = log_slope * log10(energies / energies(1)) + 1;
expected = 1 + floor((log10(test_e / energies(1))) * log_slope);

figure(1)
semilogx(energies,num, energies,y)
hold on
semilogx(test_e,expected,'o')
hold off

for ii = 1:length(energies)
    [s, l, t, S] = cs.CrossSectionCalc_Norm(Z, energies(ii), ion_level, w_p);
    
    sig(ii) = s;
    lost(ii) = l;
    transf(ii) = t;
    shell_sig(ii,:) = S;
end

% lin_s = lin_inter_extra_polate(1, 2, test_e, energies, sig);
% log_s = log_inter_extra_polate(1, 2, test_e, energies, sig);
[act_s, act_l, act_t, ~] = cs.CrossSectionCalc_Norm(Z, test_e, ion_level, w_p);

% test_e = 177;
[sigma, e_lost, e_trans] = get_value_interpolation_search(test_e, cs, Z, ion_level, occ_idx, num_points, log_slope, energies, sig, lost, transf);

% extrap_size = 50;
% extrap_energies = logspace(log10(e_max),log10(5*test_e),extrap_size);
% extrap_sig = zeros(size(extrap_energies));
% extrap_lost = zeros(size(extrap_energies));
% extrap_transf = zeros(size(extrap_energies));
% extrap_shell_sig = zeros(length(extrap_energies),occ_idx);
% for ii = 1:extrap_size
%     [s, l, t, S] = cs.CrossSectionCalc(Z, extrap_energies(ii), ion_level);
%     
%     extrap_sig(ii) = s;
%     extrap_lost(ii) = l;
%     extrap_transf(ii) = t;
%     extrap_shell_sig(ii,:) = S;
% end
% 
figure(2)
plot_title = ['\textbf{Z = ',num2str(Z),'}'];
sgtitle(plot_title, 'Interpreter','latex')
set(gca,'FontSize',20)

subplot(3,1,1)
semilogx(energies,sig, 'LineWidth',1.5)
% loglog(energies,sig, 'LineWidth',1.5)
hold on
semilogx(test_e,act_s,'o')
semilogx(test_e,sigma,'*')
% semilogx(test_e,lin_s,'+')
% semilogx(test_e,log_s,'*')
% semilogx(extrap_energies,extrap_sig, 'LineWidth',1.5)
% semilogx(test_e,sigma,'o')
% xline(test_e,'LineWidth',1.5)
% yline(sigma,'LineWidth',1.5)
% for ii = occ_idx:-1:1
%     loglog(energies,shell_sig(:,ii), 'LineWidth',1.5)
% end
% xline(0.414*mec2, 'LineWidth',1.5)
% xline(4e9, 'LineWidth',1.5)
title('\textbf{Cross Section}', 'Interpreter','latex')
xlabel('\textbf{T (eV)}', 'Interpreter','latex')
ylabel('\textbf{$\sigma$ (m$^2$)}', 'Interpreter','latex')
legend('','Actual Val','Linear Interp','Log Interp', 'location','best')
hold off
set(gca,'FontSize',16)

subplot(3,1,2)
semilogx(energies,lost, 'LineWidth',1.5)
% loglog(energies,lost)
hold on
semilogx(test_e,act_l,'o')
% semilogx(extrap_energies,extrap_lost, 'LineWidth',1.5)
% semilogx(test_e,e_lost,'o')
% xline(test_e,'LineWidth',1.5)
% yline(e_lost,'LineWidth',1.5)
title('\textbf{Energy Lost}', 'Interpreter','latex')
xlabel('\textbf{T (eV)}', 'Interpreter','latex')
ylabel('\textbf{$w$ (eV)}', 'Interpreter','latex')
hold off
set(gca,'FontSize',16)

subplot(3,1,3)
semilogx(energies,transf, 'LineWidth',1.5)
% loglog(energies,transf)
hold on
semilogx(test_e,act_t,'o')
% semilogx(extrap_energies,extrap_transf, 'LineWidth',1.5)
% semilogx(test_e,e_trans,'o')
% xline(test_e,'LineWidth',1.5)
% yline(e_trans,'LineWidth',1.5)
title('\textbf{Energy Transferred}', 'Interpreter','latex')
xlabel('\textbf{T (eV)}', 'Interpreter','latex')
ylabel('\textbf{$\varepsilon$ (eV)}', 'Interpreter','latex')
hold off
set(gca,'FontSize',16)


%% Search and Interp function
%--------------------------------------------------------------------------
% Get the cross section, energy transferred, and energy lost from table by
% interpolating from the precomputed values
%--------------------------------------------------------------------------
function [sigma, e_lost, e_trans] = get_value_interpolation_search(e_energy, cs, Z, ion_level, occ_idx, num_points, slope, energies, sig, lost, transf)
    np = num_points;

    if (e_energy < energies(1))
        if (e_energy < cs.GetBindingEnergy(Z, ion_level, occ_idx)) %? %%%(this%n_occ(ion_state)))
            sigma   = 0;
            e_lost  = 0;
            e_trans = 0;
            return
        else
            idx = 1;
        end
    elseif (e_energy > energies(np))
        idx = np - 1;
    else
        idx = interpolation_search_SIP(energies, e_energy, slope);

        if (energies(idx) == e_energy)
            % values already calculated
            sigma = sig(idx);
            e_lost = lost(idx);
            e_trans = transf(idx);
            return
        end
    end
    
    sigma = log_inter_extra_polate(idx, idx+1, e_energy, energies, sig);
    
    e_lost = log_inter_extra_polate(idx, idx+1, e_energy, energies, lost);
    
    e_trans = log_inter_extra_polate(idx, idx+1, e_energy, energies, transf);
    
    if (sigma <= 0.0 || e_lost <= 0.0 || e_trans <= 0.0)
        sigma = 0.0;
        e_lost = 0.0;
        e_trans = 0.0;
    end
end


%% Interpolation functions
%--------------------------------------------------------------------------
% Function to linear interpolate / extrapolate values
%--------------------------------------------------------------------------
function y = lin_inter_extra_polate(idx1, idx2, x_ext, x_arr, y_arr)
    a = x_ext - x_arr(idx1);
    b = x_arr(idx2) - x_ext;
    f = a / (a + b);

    y = ((1.0 - f) * y_arr(idx1)) + (f * y_arr(idx2));
end

%--------------------------------------------------------------------------
% Function to logrithmically interpolate / extrapolate values
%--------------------------------------------------------------------------
function y = log_inter_extra_polate(idx1, idx2, x_ext, x_arr, y_arr)
    a = log10(x_ext / x_arr(idx1));
    b = log10(x_arr(idx2) / x_ext);
    f = a / (a + b);

    y = y_arr(idx2).^f * y_arr(idx1).^(1.0 - f);
end


%% Search functions
%--------------------------------------------------------------------------
% SIP (slope reuse interpolation search) algorithm
% assumes a sorted array to search through, and that the value being searched
% for is within the min and max values of the array
% this should be faster than bisection search, especially since the energy
% array we search through is linearly distributed
% algorithm modified to be "fuzzy" from
% P. Van Sandt, Y. Chronis, and J. M. Patel, in Proceedings of The 2019
% International Conference on Management of Data - SIGMOD ’19 (ACM Press,
% Amsterdam, Netherlands, 2019), pp. 36–53.
%--------------------------------------------------------------------------
function idx = interpolation_search_SIP(array, to_find, slope)
    guard_size = 8;  % recommended value from paper

    n = length(array);
    left = 1;
    right = n;

    temp = floor((log10(to_find / array(left))) * slope);
%     temp = (to_find - array(left)) * slope;
    idx = left + floor(temp);  % use floor here to 'cast as int'

    % handle edge cases where the data is perfectly linear and we look at the
    % extremes of the distribution. otherwise, we deal with going out of bounds
    % on the array
    if (idx == 1 || idx == n)
        return
    end

    while ( true )
        % idx is always a lower bound since we round down going to integers
        % thus, we only have to check the cell to the right to see if we found
        % right index
        if (array(idx) <= to_find && array(idx+1) >= to_find)
            return  % idx is already the correct value
        elseif (array(idx) < to_find)
            left = idx + 1;
        elseif (array(idx) > to_find)
            right = idx - 1;
        end

        temp = floor((log10(to_find / array(idx))) * slope);
%         temp = (to_find - array(idx)) * slope;
        idx = idx + floor(temp);  % use floor here to 'cast as int'

        if (idx + guard_size >= right)
            idx = sequential_search(array, to_find, right, -1);
            return
        elseif (idx - guard_size <= left)
            idx = sequential_search(array, to_find, left, 1);
            return
        end
    end
end

%--------------------------------------------------------------------------
% Linearly search either up or down the array for a value
%--------------------------------------------------------------------------
function idx = sequential_search(array, to_find, start_idx, dir)
    idx = start_idx;

    if (dir < 0)
        while (idx >= 1)
            if (array(idx) >= to_find && array(idx-1) <= to_find)
                idx = idx - 1;  % want the index to the left
                return
            else
                idx = idx + dir;
            end
        end
    else
        while (idx <= size(array))
            if (array(idx) <= to_find && array(idx+1) >= to_find)
                return  % idx should already be the correct value
            else
                idx = idx + dir;
            end
        end
    end
end

function idx = bisection_search(array, to_find)
    n = length(array);
    % could require that the value to look for be within bounds before
    % doing search like above. left in for now in case it's used elsewhere
    if (to_find < array(1))
        idx = -1;
        return
    elseif (to_find > array(n))
        idx = n;
        return
    end

    l = 1;
    u = n;
    while (u - l > 1)
        m = bitshift(u + l, -1);
        if (to_find >= array(m))
            l = m;
        else
            u = m;
        end
    end

    if (to_find == array(1))
        idx = 1;
        return
    elseif (to_find == array(n))
        idx = n;
        return
    else
        idx = l;
        return
    end
end
