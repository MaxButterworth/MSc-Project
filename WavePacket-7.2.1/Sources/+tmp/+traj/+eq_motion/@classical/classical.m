%--------------------------------------------------------------------------
%
% Fully classical trajectories
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

classdef classical < handle
    
    properties (Access = public)
    end    

    methods (Access = public)
        
        % Construct object
        function obj = classical
        end
            
        % Initialize
        function eq_init (obj, state, propa_type)   
            init ( propa_type, state, obj)
        end
        
        % Propagate
        function eq_propa (obj, state, propa_type) 
            propa ( propa_type, state, obj)
        end
        
        % Optionally solve TDSEs attached to trajectories
        function psi = solve_tdse_adi (obj, state, step_size, mom, pot_mats, frc_mats, D, U, mom_old, pot_mats_old, frc_mats_old, D_old, U_old, psi_old, psi_oldold)
            psi = psi_old;
            if (state.extra_tdse_solving)
                switch state.tdse_solver
                    case 1
                        psi = tdse_adi (state, step_size, D, U, U_old, psi_old);
                    case 2
                        psi = tdse_adi_euler (state, step_size, mom, pot_mats, frc_mats, D, U, psi_old, psi_oldold);
                    case 3
                        psi = tdse_adi_trapezoidal (state, step_size, mom, pot_mats, frc_mats, D, U, mom_old, pot_mats_old, frc_mats_old, D_old, U_old, psi_old, psi_oldold);
                end
            end
        end
        
        % Optionally solve TDSEs attached to trajectories
        function psi = solve_tdse_dia (obj, state, step_size, pos, psi_old, psi_oldold)
            psi = psi_old;
            if (state.extra_tdse_solving)
                switch state.tdse_solver
                    case 1
                        psi = tdse_dia_2 (state, step_size, pos, psi_old, psi_oldold);
                    case 2
                        psi = tdse_dia (state, step_size, pos, psi_old, psi_oldold);
                end
            end
        end
        
        % Calculate forces in adiabatic representation
        function frc_adi = eq_frc (obj, pos, mom, pot_mats, frc_mats, D, U, m, psi) 
            frc_adi = ham.frc_adi (pot_mats,frc_mats,U,m);
        end
        
        % Calculate forces in diabatic representation
        function frc_dia_m = eq_frc_dia (obj, pos, m, psi)
            frc_dia_m = ham.frc_dia_m(pos, m);
        end
        
        
        % see separate files for the following public methods
        eval_V_F  ( obj, state, step_size, init)% Evaluate potential energy and forces 
        
    end
end



    
