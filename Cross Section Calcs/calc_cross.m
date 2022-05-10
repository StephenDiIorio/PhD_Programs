%% Cross Sectioz Set Up
% syms T

c = 299792458;
mec2 = 510998.95;
r_e = 2.8179403262E-15;
eps_0 = 8.8541878128E-12;
e_charge = 1.602176634E-19;
me = 9.1093837015E-31;
n_0 = 2.16e23 * 1e6;
w_p = sqrt(n_0 * e_charge * e_charge / me / eps_0);
toJ = 1 / 6.24e18;

Z = 8;
ion_level = 3;

cs = CrossSection;

nocc = cs.GetMaxOccShell(Z, ion_level);


%% Cross Section Test
% energies = linspace(cs.GetMinBindingEnergy(Z, ion_level), 1e6, 100000);
% energies = linspace(1e4, 1e9, 100000);
% energies = logspace(0,12,1000);
energies = logspace(log10(cs.GetMinBindingEnergy(Z, ion_level)+0.001),log10(1e12),1000);
sig = zeros(size(energies));
lost = zeros(size(energies));
transf = zeros(size(energies));
shell_sig = zeros(length(energies),nocc);
shell_TT = zeros(length(energies),nocc);

moller_sig = zeros(size(energies));

sig_norm = zeros(size(energies));
lost_norm = zeros(size(energies));
transf_norm = zeros(size(energies));
shell_sig_norm = zeros(length(energies),nocc);
shell_TT_norm = zeros(length(energies),nocc);
for ii = 1:length(energies)
    [s, l, t, S, TT] = cs.CrossSectionCalc(Z, energies(ii), ion_level);
    [sn, ln, tn, Sn, TTn] = cs.CrossSectionCalc_Norm(Z, energies(ii)/mec2, ion_level, w_p);
%     [sm, ~] = cs.MollerCrossSectionCalc(Z, energies(ii), ion_level);
    sig(ii) = s;
    lost(ii) = l;
    transf(ii) = t;
    shell_sig(ii,:) = S;

    shell_TT(ii,:) = TT;

%     moller_sig(ii) = sm;
    
    sig_norm(ii) = sn;
    lost_norm(ii) = ln;
    transf_norm(ii) = tn;
    shell_sig_norm(ii,:) = Sn;
    shell_TT_norm(ii,:) = TTn;
end
[s, ~, ~, ~] = cs.CrossSectionCalc(Z, 0.414*mec2, ion_level);
[s, ~, ~, ~] = cs.CrossSectionCalc_Norm(Z, 0.414, ion_level, w_p);
% [s, ~, ~, ~] = cs.CrossSectionCalc(Z, 4e9, ion_level);
% [s, ~, ~, ~] = cs.CrossSectionCalc(Z, 0.414*mec2, ion_level+1);
% [s, ~, ~, ~] = cs.CrossSectionCalc_Norm(Z, 0.414, ion_level+1, w_p);
% [s, ~, ~, ~] = cs.CrossSectionCalc(Z, 4e9, ion_level+1);

if ion_level-1 == 0
    plot_title = ['\textbf{Z = ',num2str(Z),', Neutral Atom}'];
else
    plot_title = ['\textbf{Z = ',num2str(Z),', Ionization Level = +',num2str(ion_level-1),'}'];
end

figure(1)
% plot_title = ['\textbf{Neutral Z = ',num2str(Z),'}'];
sgtitle(plot_title, 'Interpreter','latex')
set(gca,'FontSize',20)

subplot(3,1,1)
% semilogx(energies,sig)
loglog(energies,sig_norm, 'LineWidth',1.5)
hold on
% loglog(energies/mec2,moller_sig, 'LineWidth',1.5)
for ii = nocc:-1:1
    loglog(energies,shell_sig_norm(:,ii), 'LineWidth',1.5)
end
% xline(50e3/mec2, 'LineWidth',1.5)
% vel = 1.3173979779892139E-002;
% eng = -(-1 + vel^2 + sqrt(1 - vel^2))/(-1 + vel^2);
% xline(eng, 'LineWidth',1.5)
% vel = 1.7753230244673377E-002;
% eng = -(-1 + vel^2 + sqrt(1 - vel^2))/(-1 + vel^2);
% xline(eng, 'LineWidth',1.5)
% xline(4e9, 'LineWidth',1.5)
title('\textbf{Total Cross Section}', 'Interpreter','latex')
xlabel('\textbf{T (eV)}', 'Interpreter','latex')
ylabel('\textbf{$\sigma$ (m$^2$)}', 'Interpreter','latex')
hold off
set(gca,'FontSize',16)

subplot(3,1,2)
semilogx(energies,lost_norm*mec2, 'LineWidth',1.5)
% loglog(energies,lost)
title('\textbf{Incident Electron Energy Lost}', 'Interpreter','latex')
xlabel('\textbf{T (eV)}', 'Interpreter','latex')
ylabel('\textbf{$\varepsilon$ (eV)}', 'Interpreter','latex')
set(gca,'FontSize',16)

subplot(3,1,3)
semilogx(energies,transf_norm*mec2, 'LineWidth',1.5)
% loglog(energies,transf)
hold on
for ii = nocc:-1:1
    semilogx(energies,shell_TT_norm(:,ii)*mec2, 'LineWidth',1.5)
end
hold off
title('\textbf{Energy Transferred to Ionized Electron}', 'Interpreter','latex')
xlabel('\textbf{T (eV)}', 'Interpreter','latex')
ylabel('\textbf{$w$ (eV)}', 'Interpreter','latex')
set(gca,'FontSize',16)

% figure(3)
% subplot_title = ['\textbf{Z = ',num2str(Z),', K-Shell}'];
% semilogx(energies,shell_a1(:,1), energies,shell_a2(:,1), energies,shell_a3(:,1), energies,shell_test(:,1), 'LineWidth',1.5)
% title(subplot_title, 'Interpreter','latex')
% xlabel('\textbf{T (eV)}', 'Interpreter','latex')
% legend('$A_1$','$A_2$','$A_3$','$\ln(\frac{\beta_t^2}{(1.0 - \beta_t^2)})$', 'Location','best', 'Interpreter','latex')
% set(gca,'FontSize',16)

figure(2)
sgtitle(plot_title, 'Interpreter','latex')
set(gca,'FontSize',20)

subplot(3,1,1)
semilogx(energies, sig_norm.*(c/w_p)^2, 'LineWidth',1.5)
title('\textbf{Total Cross Section}', 'Interpreter','latex')
xlabel('\textbf{T (eV)}', 'Interpreter','latex')
ylabel('\textbf{$\sigma$ (m$^2$)}', 'Interpreter','latex')
hold off
set(gca,'FontSize',16)

subplot(3,1,2)
semilogx(energies, lost_norm.*mec2, 'LineWidth',1.5)
title('\textbf{Incident Electron Energy Lost}', 'Interpreter','latex')
xlabel('\textbf{T (eV)}', 'Interpreter','latex')
ylabel('\textbf{$\varepsilon$ (eV)}', 'Interpreter','latex')
set(gca,'FontSize',16)

subplot(3,1,3)
semilogx(energies, transf_norm.*mec2, 'LineWidth',1.5)
title('\textbf{Energy Transferred to Ionized Electron}', 'Interpreter','latex')
xlabel('\textbf{T (eV)}', 'Interpreter','latex')
ylabel('\textbf{$w$ (eV)}', 'Interpreter','latex')
set(gca,'FontSize',16)

% T = 100;
% [s, l, t, S] = CrossSectionCalc(T, ion_level, n_occ, binding_energy, xi2, e_occupation);
