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

classdef fssh_3 < sht.mssh & handle
    
    methods (Access = public)

        % Constructor: Setting defaults and a text string
        function obj = fssh_3 (n,seed)
            obj = obj@sht.mssh(n,seed);     % Inherit from superclass
            obj.string0 = 'fssh_3';
            obj.string4 = 'Fewest switches surface hopping 3';
        end   
                
        % Get probabilities of hopping from state "m" to state "n"
        function probable = prob_hop (obj,m,n,ind_m)
            global hamilt
            probable = zeros(size(ind_m));
            
            c_old_m = obj.psi_old{m}(ind_m);
            rho_m = abs(c_old_m).^2;
            rho_n = abs(obj.psi_old{n}(ind_m)).^2;
            
            nominator = zeros(size(ind_m));
            
            rho_nl = cell(hamilt.coupling.n_eqs,1);
            
            for l = 1 : hamilt.coupling.n_eqs
                rho_nl{l}  = zeros(size(ind_m));
                for k = 1 : hamilt.coupling.n_eqs
                    rho_nl{l}(:) = rho_nl{l} + squeeze( conj(obj.U_new(k,n,ind_m)) .* obj.U_old(k,l,ind_m) );
                end
                nominator(:) = nominator + abs(obj.psi_old{l}(ind_m)).^2 .* abs(rho_nl{l}).^2;
            end
            
            nominator = nominator + 2 .* real( conj(obj.psi_old{1}(ind_m)) .* obj.psi_old{2}(ind_m) .* conj(rho_nl{1}) .* rho_nl{2} );
            
            p =  ( nominator - rho_n ) ./ rho_m;
            
            probable = max( p , probable );
            
        end
        
    end
end

