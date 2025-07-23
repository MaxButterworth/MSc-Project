%--------------------------------------------------------------------------
%
% Quantum trajectories in adiabatic representation
% with non-classical contributions to forces
%
% see work by Craig C. Martens
% DOI:10.1021/acs.jpca.8b10487
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2020-.... Leonardo Cancissu Araujo
%
% see the README file for license details.
%  comment added by CCM on May 20, 2020

classdef quantraj_adi_L_mod < tmp.traj.eq_motion.classical & handle
    
    properties (Access = public)
        add_bmm
    end    

    methods (Access = public)
        
        % Construct object
        function obj = quantraj_adi_L_mod
            global hamilt
            
            if(~strcmpi(hamilt.coupling.represent,'adi'))
                prt.error('Only for trajectories following adiabatic surfaces');
            end
            
            obj.add_bmm = 0;
        end
        
        % Initialize
        function eq_init (obj, state, propa_type) 
            
%             if ~strcmp(propa_type.string,'TBD...')
%                 prt.error('qt_adi can only be used for propa = "TBD...", i.e. set qm_propa("TBD...")  ')
%             end
            
            eq_init@tmp.traj.eq_motion.classical ( obj, state, propa_type)
            
            state.extra_tdse_solving = 1;
            
        end
        
        
        % Propagate
        function eq_propa (obj, state, propa_type) 
            global space time hamilt
            
            extra_dpos_dt = cell(space.n_dim,1);
            for d=1:space.n_dim
                extra_dpos_dt{d} = zeros(length(state.pos{1}),1);
            end
            for m=1:hamilt.coupling.n_eqs
                % Indices of trajectories in that channel
                ind_m = find(state.cha==m);
                
                if(~isempty(ind_m))
                    
                    frc_mats_ind_m = cell(space.n_dim,1);
                    pos_ind_m      = cell(space.n_dim,1);
                    mom_ind_m      = cell(space.n_dim,1);
                    
                    for d=1:space.n_dim
                        frc_mats_ind_m{d}   = state.frc_mat{d}(:,:,ind_m);
                        pos_ind_m{d}        = state.pos{d}(ind_m);
                        mom_ind_m{d}        = state.mom{d}(ind_m);
                    end
                    
                    psi_ind_m = cell(hamilt.coupling.n_eqs,1);
                    
                    for n=1:hamilt.coupling.n_eqs
                        psi_ind_m{n} = state.psi{n} (ind_m);
                    end
                    
                    % Get extra term for propagation of position
                    % NAC vectors are twice: here and when getting the force
                    % -> NAC needs to become a property of traj!
                    extra_dpos_dt_m = get_extra_dpos_dt (obj, pos_ind_m, mom_ind_m, state.pot_mat(:,:,ind_m), ...
                        frc_mats_ind_m, state.D_new(:,ind_m), state.U_new(:,:,ind_m), psi_ind_m, m);
                    
                    for d=1:space.n_dim
                        extra_dpos_dt{d}(ind_m) = extra_dpos_dt_m{d};
                    end
                end
            end
            
            % Propagate positions by full time step
            for d = 1:space.n_dim
                state.pos{d} = state.pos{d} + time.steps.s_delta * state.mom{d} / space.dof{d}.mass ...
                                            + time.steps.s_delta * extra_dpos_dt{d};
            end
            
            % Get potentials and forces at new positions
            eval_V_F (obj, state, time.steps.s_delta, 0)
            
            % Propagate momenta by full time step
            for d = 1:space.n_dim
                state.mom{d} = state.mom{d} + time.steps.s_delta * state.frc{d};
            end
            
        end
        
        % Propagate
        function frc_adi = eq_frc (obj, pos, mom, pot_mats, frc_mats, D, U, m, psi) 
            global space
            
            % Get classical adiabatic force
            frc_adi = eq_frc@tmp.traj.eq_motion.classical...
                        (obj, pos, mom, pot_mats, frc_mats, D, U, m, psi);
            
            % Add non-classical contributions to adiabatic force
            extra_frc_qt_adi = get_extra_frc_qt_adi(obj, pos, mom, pot_mats, frc_mats, D, U, psi, m); 
            
            for d = 1:space.n_dim
                frc_adi{d} = frc_adi{d} + extra_frc_qt_adi{d};
            end
        end
        
        % Calculate non-classical contributions to adiabatic force
        function extra_frc_qt_adi = get_extra_frc_qt_adi(obj, pos, mom, pot_mats, frc_mats, D, U, psi, m) 
            global space hamilt
            extra_frc_qt_adi = cell(space.n_dim,1);
            
            % Compute all nac vectors - Actually, not all nac vectors are
            % needed -> performance can be improved
            nac_all = ham.nac_all (pot_mats,frc_mats,U,D);
            
            % Compute derivative of diabatic forces
            dfrc_dx_mats = ham.dfrc_dx_dia(pos);
            
            % Compute adiabatic force on each level
            frc_all = ham.frc_adi_all (pot_mats,frc_mats,U);
            
            % Compute derivative of all nac-vectors
            dnac_dx_all = ham.dnac_dx_all (dfrc_dx_mats,U,D,frc_all,nac_all);
            
            n_traj = length(pos{1});
            
            prob_all = get_fssh_prob_all(obj, pos, mom, pot_mats, frc_mats, D, U, psi, m) ;
            
            for k=1:space.n_dim
                extra_frc_k = zeros (n_traj,1);
                
                for n = 1:hamilt.coupling.n_eqs
                    if (n~=m)
                        b = imag(psi{n} ./ abs(psi{n}) .* abs(psi{m}) ./ psi{m});
                        
                        b = sqrt( prob_all{n} ./ prob_all{m} ) .* b;
                        
                        for i = 1:space.n_dim
                            extra_frc_k = extra_frc_k - b .* dnac_dx_all{m,n,k,i} .* mom{i} / space.dof{i}.mass;
                        end
                    end
                end
                extra_frc_qt_adi{k} = extra_frc_k;
            end
            
            if(obj.add_bmm)
                % Get db_mm_dx
                db_mm_dx = ham.db_mm_dx (nac_all,dnac_dx_all,m);
                for k=1:space.n_dim
                    extra_frc_qt_adi{k} = extra_frc_qt_adi{k} + db_mm_dx{k};
                end
            end
            
        end
        
        % Calculate non-classical contributions to adiabatic force
        function extra_dpos_dt = get_extra_dpos_dt(obj, pos, mom, pot_mats, frc_mats, D, U, psi, m) 
            global space hamilt
            extra_dpos_dt = cell(space.n_dim,1);
            
            % Compute all nac vectors - Actually, not all nac vectors are
            % needed -> performance can be improved
            nac_all = ham.nac_all (pot_mats,frc_mats,U,D);
            
            n_traj = length(pos{1});
            
            prob_all = get_fssh_prob_all(obj, pos, mom, pot_mats, frc_mats, D, U, psi, m);
            
            for k=1:space.n_dim
                extra_pos_k = zeros (n_traj,1);
                
                for n = 1:hamilt.coupling.n_eqs
                    if (n~=m)
                        b = imag( psi{n} ./ abs(psi{n}) .* abs(psi{m}) ./ psi{m} );
                        
                        b = sqrt( prob_all{n} ./ prob_all{m} ) .* b;
                        
                        extra_pos_k = extra_pos_k + b .* nac_all{m,n,k};
                    end
                end
                
                extra_dpos_dt{k} = extra_pos_k / space.dof{1}.mass;
            end
            
        end
        
        % Get FSSH probabilities of hopping from state "m" to state "n"
        % Eqs. (14,15) from doi:10.1063/1.459170
        function prob_mn = prob_hop_fssh_adi (obj,mom,pot_mats,frc_mats,U,D,psi,m,n)
            global space time
            
            % Preallocate
            n_traj = length(psi{1});
            prob_mn = zeros(n_traj,1);
            
            % Quantum coherence
            coherence = conj ( psi{n} ) .* psi{m};
            
            % Calculate NAC vector
            nac_mn = ham.nac_mn (pot_mats,frc_mats,U,D,m,n);
            
            % Calculate coupling (see FSSH formula)
            coupling = zeros (size(prob_mn));
            for d = 1:space.n_dim
                coupling = coupling + nac_mn{d} .* mom{d} / space.dof{d}.mass;
            end
            
            % Compute FSSH formula
            prob_mn = 2 * real (coupling .* coherence);
            
            prob_mn = time.steps.s_delta * prob_mn ./ abs(psi{m}).^2;
            
        end
        
        function prob_all = get_fssh_prob_all(obj, pos, mom, pot_mats, frc_mats, D, U, psi, m) 
            global hamilt
            n_traj = length(pos{1});
            
            prob_all = cell(hamilt.coupling.n_eqs,1);
            sum_prob = zeros(n_traj,1);
            
            for n = 1:hamilt.coupling.n_eqs
                if (n~=m)
                    prob_mn = prob_hop_fssh_adi (obj,mom,pot_mats,frc_mats,U,D,psi,m,n);
                    prob_mn = max( prob_mn, 0 );
                    prob_all{n} = prob_mn;
                    sum_prob    = sum_prob + prob_mn;
                end
            end
            
            prob_all{m} = max(1 - sum_prob,10^-5);
            
            sum_prob = sum_prob + prob_all{m};
            
            for n = 1:hamilt.coupling.n_eqs
                if 1 %(n~=m)
                    prob_all{n} = prob_all{n} ./ sum_prob;
                end
            end
            
            for n = 1:hamilt.coupling.n_eqs
                a = prob_all{n} ./ prob_all{m};
                b = a(a>1);
                if ( isempty(b) == 0 )
                    1;
                end
            end
        end
        
    end
end



    
