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

classdef quantraj_adi_L < tmp.traj.eq_motion.classical & handle
    
    properties (Access = public)
        add_bmm
    end    

    methods (Access = public)
        
        % Construct object
        function obj = quantraj_adi_L
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
            global space time
            
            % Get extra term for propagation of position
            % NAC vectors are twice calculated: here and when getting the force
            % -> NAC needs to become a property of traj!
            extra_dpos_dt = get_extra_dpos_dt (obj, state.pos, state.pot_mat, state.frc_mat,...
                                                state.D_new, state.U_new, state.psi);
            
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
        
        % Get force
        function frc_adi = eq_frc (obj, pos, mom, pot_mats, frc_mats, D, U, m, psi) 
            global space
            
            % Get classical adiabatic force
            frc_adi = eq_frc@tmp.traj.eq_motion.classical...
                        (obj, pos, mom, pot_mats, frc_mats, D, U, m, psi);
            
            % Add non-classical contributions to adiabatic force
            extra_frc_qt_adi = get_extra_frc_qt_adi(obj, pos, mom, pot_mats, frc_mats, D, U, psi); 
            
            for d = 1:space.n_dim
                frc_adi{d} = frc_adi{d} + extra_frc_qt_adi{d};
            end
        end
        
        % Calculate non-classical contributions to adiabatic force
        function extra_frc_qt_adi = get_extra_frc_qt_adi(obj, pos, mom, pot_mats, frc_mats, D, U, psi) 
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
            
            for k=1:space.n_dim
                extra_frc_k = zeros (n_traj,1);
                
                for m = 2:hamilt.coupling.n_eqs
                    for n = 1:(m-1)
                        coherence = conj ( psi{m} ) .* psi{n};
                        b         = imag(coherence);
                        
                        for i = 1:space.n_dim
                            extra_frc_k = extra_frc_k - 2.* b .* dnac_dx_all{m,n,k,i} .* mom{i} / space.dof{i}.mass;
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
        
        % Calculate non-classical contributions to the momentum
        function extra_dpos_dt = get_extra_dpos_dt(obj, pos, pot_mats, frc_mats, D, U, psi) 
            global space hamilt
            extra_dpos_dt = cell(space.n_dim,1);
            
            % Compute all nac vectors - Actually, not all nac vectors are
            % needed -> performance can be improved
            nac_all = ham.nac_all (pot_mats,frc_mats,U,D);
            
            n_traj = length(pos{1});
            
            for k=1:space.n_dim
                extra_pos_k = zeros (n_traj,1);
                
                for m = 2:hamilt.coupling.n_eqs
                    for n = 1:(m-1)
                        coherence = conj ( psi{m} ) .* psi{n};
                        b         = imag(coherence);
                        
                        extra_pos_k = extra_pos_k + 2.* b .* nac_all{m,n,k};
                    end
                end
                
                extra_dpos_dt{k} = extra_pos_k / space.dof{k}.mass;
            end
            
        end
        
    end
end



    
