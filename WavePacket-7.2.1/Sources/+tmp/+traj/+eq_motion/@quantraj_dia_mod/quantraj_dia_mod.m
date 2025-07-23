%--------------------------------------------------------------------------
%
% Quantum trajectories in diabatic representation
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

classdef quantraj_dia_mod < tmp.traj.eq_motion.classical & handle
    
    properties (Access = public)
    end    

    methods (Access = public)
        
        % Construct object
        function obj = quantraj_dia_mod
            global hamilt
            
            if(~strcmpi(hamilt.coupling.represent,'dia'))
                prt.error('Only for trajectories following diabatic surfaces');
            end
        end
        
        % Initialize
        function eq_init (obj, state, propa_type) 
            
            eq_init@tmp.traj.eq_motion.classical ( obj, state, propa_type)
            
            state.extra_tdse_solving = 1;
            
        end
        
        % Diabatic quantraj force
        function frc_dia_m = eq_frc_dia (obj, pos, m, psi)
            global space 
            
            % Get classical diabatic force
            frc_dia_m = ham.frc_dia_m(pos, m);
            
            % Add non-classical contributions to diabatic force
            extra_frc_qt_dia = get_extra_frc_qt_dia(obj, pos, psi, m); 
            
            for d = 1:space.n_dim
                frc_dia_m{d} = frc_dia_m{d} + extra_frc_qt_dia{d};
            end
            
        end
        
        % Extra term of diabtic quantraj force
        function extra_frc_qt_dia = get_extra_frc_qt_dia (obj, pos, psi, m)
            global hamilt space
            extra_frc_qt_dia = cell(space.n_dim,1);
            n_traj = length(pos{1});
            
            frc_mats = ham.frc_dia(pos);
            
            pot_mats = ham.pot_dia(pos);
%             prob_all = get_fssh_prob_all(obj, pos, pot_mats, psi, m);
            
            for k=1:space.n_dim
                extra_frc_k = zeros (n_traj,1);
                
                for n = 1:hamilt.coupling.n_eqs
                    if (n~=m)
                        a = real( psi{n} ./ psi{m} );
                        
%                         a = sqrt( prob_all{n} ./ prob_all{m} );
                        
                        frc_mn_k    = zeros (n_traj,1);
                        frc_mn_k(:) = frc_mats{k}(m,n,:);
                        
                        extra_frc_k = extra_frc_k + a .* frc_mn_k;
                        
                    end
                end
                extra_frc_qt_dia{k} = extra_frc_k;
            end
            
        end
        
        % Get FSSH probabilities of hopping from state "m" to state "n"
        % Eqs. (14,15) from doi:10.1063/1.459170
        function prob_mn = prob_hop_fssh_dia (obj,pot_mats,psi,m,n)
            global time
            
            % Preallocate
            n_traj = length(psi{1});
            prob_mn = zeros(n_traj,1);
            coupling = zeros(n_traj,1);
            
            % Quantum coherence
            coherence = conj ( psi{n} ) .* psi{m};
            
            coupling(:) = pot_mats(m,n,:);
            prob_mn(:) = +2 * imag (coupling .* coherence);
            
            prob_mn = time.steps.s_delta * prob_mn ./ abs(psi{m}).^2;
            
        end
        
        function prob_all = get_fssh_prob_all(obj, pos, pot_mats, psi, m) 
            global hamilt
            n_traj = length(pos{1});
            
            prob_all = cell(hamilt.coupling.n_eqs,1);
            sum_prob = zeros(n_traj,1);
            
            for n = 1:hamilt.coupling.n_eqs
                if (n~=m)
                    prob_mn = prob_hop_fssh_dia (obj,pot_mats,psi,m,n);
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



    
