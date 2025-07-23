%--------------------------------------------------------------------------
%
% Surface hopping with cumulative probabilities
% doi.org/10.1063/5.0024372
% doi.org/10.1063/1.459170
% See Ticket #259  
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef fssh_i < sht.mssh & handle
    
    methods (Access = public)

        % Constructor: Setting defaults and a text string
        function obj = fssh_i (n,seed)
            obj = obj@sht.mssh(n,seed);     % Inherit from superclass
            obj.string0 = 'fssh_i';
            obj.string4 = 'Fewest switches surface hopping instantaneous';
        end   
                
        % Get probabilities of hopping from state "m" to state "n"
        % Eqs. (14,15) from doi:10.1063/1.459170
        function probable = prob_hop (obj,m,n,ind_m)
            global time
            
            rho_m_old   = abs(obj.psi_old{m}(ind_m)).^2;
            rho_n_old   = abs(obj.psi_old{n}(ind_m)).^2;
            rho_n       = abs(obj.psi    {n}(ind_m)).^2;
                
            p = ( rho_n - rho_n_old ) ./ rho_m_old;
            
            p = p .* s(obj,time.steps.s_delta); % only implemented for 2 states
            
            probable = max( p , zeros(size(ind_m)) );
            
        end
        
        function prob = s(~,x)
            
            if(abs(x) > 10^(-3))
                % Use correct function if x far from 0
                prob = - expm1(-x) / x;
            else
                % Due to numerical issues, use Taylor approx close to 0
                prob = 0;
                for k=0:4
                    prob = prob + (-x)^k / factorial(k+1);
                end
            end
        end
        
    end
end

