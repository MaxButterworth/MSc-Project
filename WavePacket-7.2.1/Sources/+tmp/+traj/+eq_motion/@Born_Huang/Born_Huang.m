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

classdef Born_Huang < tmp.traj.eq_motion.classical & handle
    
    properties (Access = public)
    end    

    methods (Access = public)
        
        % Propagate
        function frc_BH = eq_frc (obj, pos, mom, pot_mats, frc_mats, D, U, m, psi) 
            global space

            % Get classical adiabatic force
            frc_BH = eq_frc@tmp.traj.eq_motion.classical...
                        (obj, pos, mom, pot_mats, frc_mats, D, U, m, psi);
            
            % Add second order non-adiabatic coupling
            frc_extra = Born_Huang_extra_frc(obj, pos, pot_mats, frc_mats, D, U, m); 
            
            for d = 1:space.n_dim
                frc_BH{d} = frc_BH{d} + frc_extra{d};
            end
        end
        
        % Calculate non-classical contributions to adiabatic force
        function frc_extra = Born_Huang_extra_frc(obj, pos, pot_mats, frc_mats, D, U, m) 
            
            % Compute all nac vectors - Actually, not all nac vectors are
            % needed -> performance can be improved
            nac_all = ham.nac_all (pot_mats,frc_mats,U,D);
            
            % Compute derivative of diabatic forces
            dfrc_dx_mats = ham.dfrc_dx_dia(pos);
            
            % Compute adiabatic force on each level
            frc_all = ham.frc_adi_all (pot_mats,frc_mats,U);
            
            % Compute derivative of all nac-vectors
            dnac_dx_all = ham.dnac_dx_all (dfrc_dx_mats,U,D,frc_all,nac_all);
            
            % Get extra force
            frc_extra = ham.db_mm_dx (nac_all,dnac_dx_all,m);
            
        end
        
    end
end



    
