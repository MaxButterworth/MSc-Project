%--------------------------------------------------------------------------
%
% FSSH = Fewest switching surface hopping algorithm
% 
% Determine transition probability from d/dt |c(t)|^2
% from TDSEs attached to each of the trajectories.
% According to the proof given by Tully, this algorithm
% really results in the lowest number of switches.
% 
% see: J. C. Tully
%      J. Chem. Phys. 93(2), 1061-1071 (1990)
%      DOI:10.1063/1.459170
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef fssh_2 < sht.mssh & handle
    
    methods (Access = public)

        % Constructor: Setting defaults and a text string
        function obj = fssh_2 (n,seed)
            obj = obj@sht.mssh(n,seed);     % Inherit from superclass
            obj.string0 = 'fssh_2';
            obj.string4 = 'Fewest switches surface hopping 2';
            obj.ens_cont = 0;               % No continuity needed
        end   
                
        % Get probabilities of hopping from state "m" to state "n"
        % Eqs. (14,15) from doi:10.1063/1.459170
        function probable = prob_hop (obj,m,n,ind_m)
            
            rho_m_old   = abs(obj.psi_old{m}(ind_m)).^2;
            rho_n_old   = abs(obj.psi_old{n}(ind_m)).^2;
            rho_n       = abs(obj.psi    {n}(ind_m)).^2;
            rho_m       = abs(obj.psi    {m}(ind_m)).^2;

            p_out = ( rho_m_old - rho_m ) ./ rho_m_old;

            p = ( rho_n - rho_n_old ) ./ rho_m_old;

            p = min( p , p_out );
            
            probable = max( p , zeros(size(ind_m)) );
            
        end
        
    end
end

