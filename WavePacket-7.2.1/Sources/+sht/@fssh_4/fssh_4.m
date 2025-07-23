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

classdef fssh_4 < sht.mssh & handle
    
    methods (Access = public)

        % Constructor: Setting defaults and a text string
        function obj = fssh_4 (n,seed)
            obj = obj@sht.mssh(n,seed);     % Inherit from superclass
            obj.string0 = 'fssh_4';
            obj.string4 = 'Fewest switches surface hopping 4';
        end   
                
        % Get probabilities of hopping from state "m" to state "n"
        function probable = prob_hop (obj,m,n,ind_m)
            global hamilt
            probable = zeros(size(ind_m));
            
            c_old_m = obj.psi_old{m}(ind_m);
            c_old_n = obj.psi_old{n}(ind_m);
            
            um_un  = zeros(size(ind_m));
            un_un  = zeros(size(ind_m));
            for k = 1 : hamilt.coupling.n_eqs
                um_un(:) = um_un + squeeze( conj(obj.U_old(k,m,ind_m)) .* obj.U_new(k,n,ind_m) );
                un_un(:) = un_un + squeeze( conj(obj.U_old(k,n,ind_m)) .* obj.U_new(k,n,ind_m) );
            end
            
            p = 2 .* real( c_old_n ./ c_old_m ) .* um_un .* un_un;
            
            probable = max( p , probable );
            
        end
        
    end
end

