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

classdef quantraj_adi_corr < tmp.traj.eq_motion.classical & handle
    
    properties (Access = public)
    end    

    methods (Access = public)
        
        % Construct object
        function obj = quantraj_adi_corr
            global hamilt
            
            if(~strcmpi(hamilt.coupling.represent,'adi'))
                prt.error('Only for trajectories following adiabatic surfaces');
            end
        end
        
        % Initialize
        function eq_init (obj, state, propa_type) 
            
            eq_init@tmp.traj.eq_motion.classical ( obj, state, propa_type)
            
            state.extra_tdse_solving = 1;
        end
        
        % Adiabatic quantraj force
        function frc_adi = eq_frc (obj, pos, mom, pot_mats, frc_mats, D, U, m, psi)
            global space 
            
            % Get classical diabatic force
            frc_adi = eq_frc@tmp.traj.eq_motion.classical...
                        (obj, pos, mom, pot_mats, frc_mats, D, U, m, psi);
            
            % Add non-classical contributions to diabatic force
            extra_frc_qt_adi = get_extra_frc_qt_adi(obj, pos, mom, pot_mats, frc_mats, D, U, psi);
            
            for d = 1:space.n_dim
                frc_adi{d} = frc_adi{d} + extra_frc_qt_adi{d};
            end
            
        end
        
        % Extra term of adiabtic quantraj force
        function extra_frc_qt_adi = get_extra_frc_qt_adi (obj, pos, mom, pot_mats, frc_mats, D, U, psi)
            global hamilt space
            extra_frc_qt_adi = cell(space.n_dim,1);
            n_traj = length(pos{1});
            
            nac_norm_all = ham.nac_norm_all(pot_mats,frc_mats,U,D);
            
            for k=1:space.n_dim
                extra_frc_k = zeros (n_traj,1);
                
                for m = 2:hamilt.coupling.n_eqs
                    for n = 1:(m-1)
                        coherence = conj ( psi{m} ) .* psi{n};
                        a         = real(coherence);
                        
                        frc_mn_k    = zeros (n_traj,1);
                        frc_mn_k(:) = nac_norm_all{m,n,k}(:);
                        
                        extra_frc_k = extra_frc_k + 2.* a .* frc_mn_k;
                        
                    end
                end
                extra_frc_qt_adi{k} = extra_frc_k;
            end
            
        end
        
    end
end



    
