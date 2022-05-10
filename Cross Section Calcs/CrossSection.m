%% Cross Section Function
classdef CrossSection
    properties
        c
        a0
        alpha
        m_ec2
        m_ec2_inv
        Rydberg_ev
        Rydberg_norm

        mrbeb_const

        bindingenergy
        bindingenergy_norm
        screeningconstants
        electron_shell_numbers
    end
    methods
        function obj = CrossSection 
            obj.c = 299792458;
            obj.a0 = 5.29177210903e-11;
            obj.alpha = 7.2973525693e-3;
            obj.m_ec2 = 0.51099895e6;
            obj.m_ec2_inv = 1.0 / obj.m_ec2;
            obj.Rydberg_ev = 13.605693122994;
            obj.Rydberg_norm = 13.605693122994 * obj.m_ec2_inv;

            obj.mrbeb_const = 2.0 * pi * obj.alpha.^4;

            load('bindingenergy.mat', 'bindingenergy');
            obj.bindingenergy = bindingenergy;
            obj.bindingenergy_norm = bindingenergy .* obj.m_ec2_inv;
            load('screeningconstants.mat', 'guerrascreeningconstants');
            obj.screeningconstants = guerrascreeningconstants;
            load('electronoccupation.mat', 'electron_shell_numbers');
            obj.electron_shell_numbers = electron_shell_numbers;
        end

        function [sigma, en_L, en_T, S, T] = CrossSectionCalc(obj, Z, e_eng, ion_level)
            [~, occ_idx] = ElectronInformation(obj, Z, ion_level);
            binding_energy = obj.bindingenergy(Z,:);

            %% Actually do calculation

            S = zeros(1, occ_idx);
            T = zeros(1, occ_idx);

            sigma = 0.0;  % [m^2]
            en_L  = 0.0;  % [eV]
            en_T  = 0.0;  % [eV]

            % make sure e_eng is at least greater than lowest binding energy
            if (e_eng >= binding_energy(occ_idx))
                for lvl_idx = 1:occ_idx
                    [sigma_k, en_lost_k, en_transfer_k] = SingleCrossSecCalc(obj, e_eng, Z, ion_level, lvl_idx);

                    sigma = sigma + sigma_k;
                    en_L = en_L + en_lost_k;
                    en_T = en_T + en_transfer_k;

                    S(lvl_idx) = sigma;
                    if (sigma > 0.0)
                        T(lvl_idx) = en_T / sigma;
                    end
                end
            end

            if (sigma > 0.0)
                en_L = en_L / sigma;
                en_T = en_T / sigma;
            end
        end
        function [sigma, en_L, en_T, S, T] = CrossSectionCalc_Norm(obj, Z, e_eng, ion_level, w_p)
            [~, occ_idx] = ElectronInformation(obj, Z, ion_level);
            binding_energy = obj.bindingenergy_norm(Z,:);

            %% Actually do calculation

            S = zeros(1, occ_idx);
            T = zeros(1, occ_idx);

            sigma = 0.0;
            en_L  = 0.0;
            en_T  = 0.0;

            % make sure e_eng is at least greater than lowest binding energy
            if (e_eng >= binding_energy(occ_idx))
                for lvl_idx = 1:occ_idx
                    [sigma_k, en_lost_k, en_transfer_k] = SingleCrossSecCalc_Norm(obj, e_eng, Z, ion_level, lvl_idx, w_p);

                    sigma = sigma + sigma_k;
                    en_L = en_L + en_lost_k;
                    en_T = en_T + en_transfer_k;

                    S(lvl_idx) = sigma;
                    if (sigma > 0.0)
                        T(lvl_idx) = en_T / sigma;
                    end
                end
            end

            if (sigma > 0.0)
                en_L = en_L / sigma;
                en_T = en_T / sigma;
            end
        end

        function [sigma_k, en_lost_k, en_transfer_k] = SingleCrossSecCalc(obj, e_eng, Z, ion_level, shell_num)
            %% WARNING: DO NOT USE en_lost_k, en_transfer_k DIRECTLY FROM THIS FUNCTION.
            %  MISSING THE 1/SIGMA TERM IN THE SUMMATION OVER ALL ELECTRON
            %  SHELLS

            %%
            [e_occupation, occ_idx] = ElectronInformation(obj, Z, ion_level);
            xi2 = obj.screeningconstants(Z,:);
            binding_energy = obj.bindingenergy(Z,:);

            if shell_num > occ_idx
                error('Atom does not have electrons in that shell');
            end

            sigma_k       = 0.0;
            en_lost_k     = 0.0;
            en_transfer_k = 0.0;

            %% Actually do calculation
            if (e_eng >= binding_energy(occ_idx))
                tp  = e_eng * obj.m_ec2_inv;
                bt2 = 1.0 - (1.0 / (1.0 + tp).^2);
                b = binding_energy(shell_num);

                if (e_eng >= b)  % nonzero only when eng is >= binding eng
                    b_inv = 1.0 / b;
                    t = e_eng * b_inv;
                    bp = b * obj.m_ec2_inv;
                    bb2 = 1.0 - (1.0 / (1.0 + bp).^2);

                    if (shell_num < occ_idx)
                        mrbeb_c = (0.3 * xi2(shell_num)) + (0.7 * xi2(shell_num + 1));
                    else  % no electrons in higher shells
                        mrbeb_c = xi2(shell_num);
                    end
                    mrbeb_c = mrbeb_c * b_inv * obj.Rydberg_ev;  % 2's in num and den cancel

                    A1 = 1.0 / (t + 1.0) * (1.0 + 2.0 * tp) / (1.0 + 0.5 * tp).^2;
                    A2 = bp * bp / (1.0 + 0.5 * tp).^2 * (t - 1.0) * 0.5;
                    A3 = log(bt2 / (1.0 - bt2)) - bt2 - log(2.0 * bp);

                    sigma_k0 = obj.mrbeb_const * obj.a0^2 * ...
                        e_occupation(ion_level, shell_num) / ...
                        (bp * (bt2 + mrbeb_c * bb2));

                    % calculate sigma, eng lost, and transferred eng for this orbital
                    % (the k terms)
                    % mrbeb cross section
                    sigma_k = sigma_k0 ...
                        * ( 0.5 * A3 ...
                        * ( 1.0 - (1.0 / (t * t)) ) + 1.0 ...
                        - (1.0 / t) + A2 - (A1 * log(t)) );
                    % pérez energy loss/transfer
                    w_k = sigma_k0 ...
                        * ( 0.5 * A3 ...
                        * ((t - 1.0).^2 / (t * (t + 1.0))) ...
                        + 2.0 * log((t + 1) * 0.5) - log(t) ...
                        + A2 * 0.25 * (t - 1.0) ...
                        - A1 * (t * log(t) - (t + 1.0) ...
                        * log(0.5 * (t + 1.0))) );

                    en_lost_k = (w_k + sigma_k) * b;
                    en_transfer_k = w_k * b;
                end
            end
        end
        function [sigma_k, en_lost_k, en_transfer_k] = SingleCrossSecCalc_Norm(obj, e_eng, Z, ion_level, shell_num, w_p)
            %% WARNING: DO NOT USE en_lost_k, en_transfer_k DIRECTLY FROM THIS FUNCTION.
            %  MISSING THE 1/SIGMA TERM IN THE SUMMATION OVER ALL ELECTRON
            %  SHELLS

            %%
            [e_occupation, occ_idx] = ElectronInformation(obj, Z, ion_level);
            xi2 = obj.screeningconstants(Z,:);
            binding_energy = obj.bindingenergy_norm(Z,:);

            if shell_num > occ_idx
                error('Atom does not have electrons in that shell');
            end

            sigma_k       = 0.0;
            en_lost_k     = 0.0;
            en_transfer_k = 0.0;

            %% Actually do calculation
            if (e_eng >= binding_energy(occ_idx))
                tp  = e_eng; % already normalized
                bt2 = 1.0 - (1.0 / (1.0 + tp).^2);
                b = binding_energy(shell_num);

                if (e_eng >= b)  % nonzero only when eng is >= binding eng
                    b_inv = 1.0 / b;
                    t = e_eng * b_inv;
                    bp = b; % already normalized
                    bb2 = 1.0 - (1.0 / (1.0 + bp).^2);

                    if (shell_num < occ_idx)
                        mrbeb_c = (0.3 * xi2(shell_num)) + (0.7 * xi2(shell_num + 1));
                    else  % no electrons in higher shells
                        mrbeb_c = xi2(shell_num);
                    end
                    mrbeb_c = mrbeb_c * b_inv * obj.Rydberg_norm;  % 2's in num and den cancel

                    A1 = 1.0 / (t + 1.0) * (1.0 + 2.0 * tp) / (1.0 + 0.5 * tp).^2;
                    A2 = bp * bp / (1.0 + 0.5 * tp).^2 * (t - 1.0) * 0.5;
                    A3 = log(bt2 / (1.0 - bt2)) - bt2 - log(2.0 * bp);

                    sigma_k0 = obj.mrbeb_const * obj.a0^2 * (w_p / obj.c)^2 * ... % normalize a0
                        e_occupation(ion_level, shell_num) / ...
                        (bp * (bt2 + mrbeb_c * bb2));

                    % calculate sigma, eng lost, and transferred eng for this orbital
                    % (the k terms)
                    % mrbeb cross section
                    sigma_k = sigma_k0 ...
                        * ( 0.5 * A3 ...
                        * ( 1.0 - (1.0 / (t * t)) ) + 1.0 ...
                        - (1.0 / t) + A2 - (A1 * log(t)) );
                    % pérez energy loss/transfer
                    w_k = sigma_k0 ...
                        * ( 0.5 * A3 ...
                        * ((t - 1.0).^2 / (t * (t + 1.0))) ...
                        + 2.0 * log((t + 1) * 0.5) - log(t) ...
                        + A2 * 0.25 * (t - 1.0) ...
                        - A1 * (t * log(t) - (t + 1.0) ...
                        * log(0.5 * (t + 1.0))) );

                    en_lost_k = (w_k + sigma_k) * b;
                    en_transfer_k = w_k * b;
                end
            end
        end

        function [sigma, S] = MollerCrossSectionCalc(obj, Z, e_eng, ion_level)
            [~, occ_idx] = ElectronInformation(obj, Z, ion_level);
            binding_energy = obj.bindingenergy(Z,:);

            %% Actually do calculation

            S = zeros(1, occ_idx);

            sigma = 0.0;  % [m^2]

            % make sure e_eng is at least greater than lowest binding energy
            if (e_eng >= binding_energy(occ_idx))
                for lvl_idx = 1:occ_idx
                    [sigma_k] = SingleMollerCrossSecCalc(obj, e_eng, Z, ion_level, lvl_idx);

                    sigma = sigma + sigma_k;

                    S(lvl_idx) = sigma;
                end
            end
        end
        function [sigma_k] = SingleMollerCrossSecCalc(obj, e_eng, Z, ion_level, shell_num)
            %% WARNING: DO NOT USE en_lost_k, en_transfer_k DIRECTLY FROM THIS FUNCTION.
            %  MISSING THE 1/SIGMA TERM IN THE SUMMATION OVER ALL ELECTRON
            %  SHELLS

            %%
            [e_occupation, occ_idx] = ElectronInformation(obj, Z, ion_level);
            binding_energy = obj.bindingenergy(Z,:);

            if shell_num > occ_idx
                error('Atom does not have electrons in that shell');
            end

            sigma_k = 0.0;

            %% Actually do calculation
            if (e_eng >= binding_energy(occ_idx))
                tp  = e_eng * obj.m_ec2_inv;
                bt2 = 1.0 - (1.0 / (1.0 + tp).^2);
                b = binding_energy(shell_num);

                if (e_eng >= b)  % nonzero only when eng is >= binding eng
                    b_inv = 1.0 / b;
                    t = e_eng * b_inv;
                    bp = b * obj.m_ec2_inv;

                    const = 4 * pi * obj.a0^2 * obj.alpha^2 ...
                        * e_occupation(ion_level, shell_num) ...
                        * (obj.Rydberg_ev * b_inv ) / bt2;

                    sigma_k = const * (1 - (1 / t) - (log(t) / (t + 1)) ...
                        * ((1 + 2 * tp) / (1 + tp)^2) + (bp^2 / (1 + tp)^2) ...
                        * ((t - 1) / 2));
                end
            end
        end


        function occ_idx = GetMaxOccShell(obj, Z, ion_level)
            [~, occ_idx] = ElectronInformation(obj, Z, ion_level);
        end

        function e = GetMinBindingEnergy(obj, Z, ion_level)
            [~, occ_idx] = ElectronInformation(obj, Z, ion_level);

            e = obj.bindingenergy(Z,occ_idx);
        end
        function e = GetMinBindingEnergy_Norm(obj, Z, ion_level)
            [~, occ_idx] = ElectronInformation(obj, Z, ion_level);

            e = obj.bindingenergy_norm(Z,occ_idx);
        end

        function e = GetBindingEnergy(obj, Z, ion_level, idx)
            [~, occ_idx] = ElectronInformation(obj, Z, ion_level);

            if idx > occ_idx
                error('No binding energy at specified shell');
            end

            e = obj.bindingenergy(Z,idx);
        end
        function e = GetBindingEnergy_Norm(obj, Z, ion_level, idx)
            [~, occ_idx] = ElectronInformation(obj, Z, ion_level);

            if idx > occ_idx
                error('No binding energy at specified shell');
            end

            e = obj.bindingenergy_norm(Z,idx);
        end

        function eqn = GetEquation(obj, Z, ion_level, shell_num)
            syms T

            %% Setup for producing formula
            [e_occupation, occ_idx] = ElectronInformation(obj, Z, ion_level);
            xi2 = obj.screeningconstants(Z,:);
            binding_energy = obj.bindingenergy(Z,:);

            if shell_num > occ_idx
                error('Atom does not have electrons in that shell');
            end

            %% Create formula
            tp  = T * obj.m_ec2_inv;
            bt2 = 1.0 - (1.0 / (1.0 + tp).^2);
            b = binding_energy(shell_num);

            b_inv = 1.0 / b;
            t = T * b_inv;
            bp = b * obj.m_ec2_inv;
            bb2 = 1.0 - (1.0 / (1.0 + bp).^2);

            if (shell_num < occ_idx)
                mrbeb_c = (0.3 * xi2(shell_num)) + (0.7 * xi2(shell_num + 1));
            else  % no electrons in higher shells
                mrbeb_c = xi2(shell_num);
            end
            mrbeb_c = mrbeb_c * b_inv * obj.Rydberg_ev;  % 2's in num and den cancel

            A1 = 1.0 / (t + 1.0) * (1.0 + 2.0 * tp) / (1.0 + 0.5 * tp).^2;
            A2 = bp * bp / (1.0 + 0.5 * tp).^2 * (t - 1.0) * 0.5;
            A3 = log(bt2 / (1.0 - bt2)) - bt2 - log(2.0 * bp);

            sigma_k0 = obj.mrbeb_const * obj.a0^2 * ...
                e_occupation(ion_level, shell_num) / ...
                (bp * (bt2 + mrbeb_c * bb2));

            % calculate sigma, eng lost, and transferred eng for this orbital
            % (the k terms)
            % mrbeb cross section
            eqn = sigma_k0 ...
                * ( 0.5 * A3 ...
                * ( 1.0 - (1.0 / (t * t)) ) + 1.0 ...
                - (1.0 / t) + A2 - (A1 * log(t)) );
        end
        function eqn = GetEquation_Norm(obj, Z, ion_level, shell_num, w_p)
            syms T

            %% Setup for producing formula
            [e_occupation, occ_idx] = ElectronInformation(obj, Z, ion_level);
            xi2 = obj.screeningconstants(Z,:);
            binding_energy = obj.bindingenergy_norm(Z,:);

            if shell_num > occ_idx
                error('Atom does not have electrons in that shell');
            end

            %% Create formula
            tp  = T; % already normalized
            bt2 = 1.0 - (1.0 / (1.0 + tp).^2);
            b = binding_energy(shell_num);

            b_inv = 1.0 / b;
            t = T * b_inv;
            bp = b; % already normalized
            bb2 = 1.0 - (1.0 / (1.0 + bp).^2);

            if (shell_num < occ_idx)
                mrbeb_c = (0.3 * xi2(shell_num)) + (0.7 * xi2(shell_num + 1));
            else  % no electrons in higher shells
                mrbeb_c = xi2(shell_num);
            end
            mrbeb_c = mrbeb_c * b_inv * obj.Rydberg_norm;  % 2's in num and den cancel

            A1 = 1.0 / (t + 1.0) * (1.0 + 2.0 * tp) / (1.0 + 0.5 * tp).^2;
            A2 = bp * bp / (1.0 + 0.5 * tp).^2 * (t - 1.0) * 0.5;
            A3 = log(bt2 / (1.0 - bt2)) - bt2 - log(2.0 * bp);

            sigma_k0 = obj.mrbeb_const * obj.a0^2 * (w_p / obj.c)^2 * ... % normalize a0
                e_occupation(ion_level, shell_num) / ...
                (bp * (bt2 + mrbeb_c * bb2));

            % calculate sigma, eng lost, and transferred eng for this orbital
            % (the k terms)
            % mrbeb cross section
            eqn = sigma_k0 ...
                * ( 0.5 * A3 ...
                * ( 1.0 - (1.0 / (t * t)) ) + 1.0 ...
                - (1.0 / t) + A2 - (A1 * log(t)) );
        end
    end


    methods(Access=private)
        function [e_occupation, occ_idx] = ElectronInformation(obj, Z, ion_level)
            e_shell_idx = 1;
            for ii = 1:Z-1
                e_shell_idx = e_shell_idx + ii;
            end

            e_occupation = obj.electron_shell_numbers(e_shell_idx:e_shell_idx+(Z-1),:);

            n_occ = zeros(1,Z);
            for ii = 1:Z  % for each ionization state
                occ_idx = 1;
                for jj = 1:size(e_occupation,2)  % for each shell
                    if (e_occupation(ii, jj) > 0)
                        occ_idx = jj;
                    end
                end
                n_occ(ii) = occ_idx;
            end

            occ_idx = n_occ(ion_level);
        end
    end
end
