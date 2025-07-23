%--------------------------------------------------------------------------
%
% Superadiabatic equations of motion
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

classdef super_adi < tmp.traj.eq_motion.classical & handle
    
    properties (Access = public)
    end    

    methods (Access = public)
        
        % Construct object
        function obj = super_adi
        end
          
        function eval_V_F (obj, state, step_size, init)
            global hamilt space time
            
            % Calculate and save all diabatic potentials in advance
            pot_mats_old = state.pot_mat;
            state.pot_mat = ham.pot_dia(state.pos);
            
            % Compute adiabatic potential matrix and eigenvector matrix
            U_old = state.U_new;
            D_old = state.D_new;
            [state.U_new , state.D_new] = ham.pot_eig_adi(state.pot_mat);
            
            % Ensure continuity of eigenvector matrix
            if( hamilt.coupling.n_eqs > 2 && state.ens_cont )
                ensure_continuity_eig (state, U_old);
            end
            
            % Calculate and save all diabatic forces in advance
            frc_mats_old = state.frc_mat;
            state.frc_mat = ham.frc_dia(state.pos);
            
            if(~init && state.extra_tdse_solving)
                psi = solve_tdse_adi (obj, state, step_size, state.mom, state.pot_mat, state.frc_mat, ...
                    state.D_new, state.U_new, ...
                    state.mom_old, pot_mats_old, frc_mats_old, D_old, U_old,...
                    state.psi, state.psi_old);

                state.psi_old = state.psi;
                state.psi     = psi;
            end
            
            % Compute all superadiabatic eigenvectors
            state.u_sa = ham.pot_eig_sa(state.mom, state.pot_mat, state.frc_mat, state.D_new, state.U_new);
            
            % Choose only the sa-eigenvectors which corresponds to the current state of the trajectory
            % u_sa_cha = state.u_sa(:,state.cha,:); 
            u_sa_cha = zeros(hamilt.coupling.n_eqs , state.n_p);
            
            % Loop over all states
            for m = 1:hamilt.coupling.n_eqs
                ind_m = find( state.cha==m );
                u_sa_cha(:,ind_m) = state.u_sa(:,m,ind_m);
            end
            
            % Reset pot and frc
            state.pot = zeros (state.n_p, 1);
            for d=1:space.n_dim
                state.frc{d} = zeros (state.n_p, 1);
            end
            
            % Loop over all states
            for m = 1:hamilt.coupling.n_eqs
                    
                % Save adiabatic potential energies
                state.pot = state.pot + abs(u_sa_cha(m,:)').^2 .* state.D_new(m,:)';
                
                frc_adi_m = ham.frc_adi (state.pot_mat,state.frc_mat,state.U_new,m);
                
                for d=1:space.n_dim
                    state.frc{d} = state.frc{d} + abs(u_sa_cha(m,:)').^2 .* frc_adi_m{d};
                end
            end
            
            time.counter.diag.dynamics = time.counter.diag.dynamics + state.n_p;
        end
        
    end
end



    
