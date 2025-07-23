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

classdef quantraj_dia < tmp.traj.eq_motion.classical & handle
    
    properties (Access = public)
    end    

    methods (Access = public)
        
        % Construct object
        function obj = quantraj_dia
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
            extra_frc_qt_dia = get_extra_frc_qt_dia(obj, pos, psi); 
            
            for d = 1:space.n_dim
                frc_dia_m{d} = frc_dia_m{d} + extra_frc_qt_dia{d};
            end
            
        end
        
        % Extra term of diabtic quantraj force
        function extra_frc_qt_dia = get_extra_frc_qt_dia (obj, pos, psi)
            global hamilt space
            extra_frc_qt_dia = cell(space.n_dim,1);
            n_traj = length(pos{1});
            
            frc_mats = ham.frc_dia(pos);
            
            for k=1:space.n_dim
                extra_frc_k = zeros (n_traj,1);
                
                for m = 2:hamilt.coupling.n_eqs
                    for n = 1:(m-1)
                        coherence = conj ( psi{m} ) .* psi{n};
                        a         = real(coherence);
                        
                        frc_mn_k    = zeros (n_traj,1);
                        frc_mn_k(:) = frc_mats{k}(m,n,:);
                        
                        extra_frc_k = extra_frc_k + 2.* a .* frc_mn_k;
                        
                    end
                end
                extra_frc_qt_dia{k} = extra_frc_k;
            end
            
        end
        
    end
end



    
