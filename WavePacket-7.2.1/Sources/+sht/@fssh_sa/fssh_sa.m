%--------------------------------------------------------------------------
% FSSH with hopping probabilities in the super adiabatic representation
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef fssh_sa < sht.mssh & handle
    
    methods (Access = public)

        % Constructor: Setting defaults and a text string
        function obj = fssh_sa (n,seed)
            obj = obj@sht.mssh(n,seed);     % Inherit from superclass
            
            obj.choice_eq_motion = 'super_adi';
            obj.rescale  = false;               % do not rescale momentum after a hop
            
            obj.string0 = 'fssh_sa';
            obj.string4 = 'Fewest switches surface hopping Super-Adi';
        end   
                
        % Get probabilities of hopping from state "m" to state "n"
        function probable = prob_hop (obj,m,n,ind_m)
            global hamilt
            probable = zeros(size(ind_m));
            
            if ~isempty(obj.u_sa_old)
                
                c_adi_old = cell(hamilt.coupling.n_eqs,1);
                c_adi     = cell(hamilt.coupling.n_eqs,1);
                for k = 1 : hamilt.coupling.n_eqs
                    c_adi_old{k} = obj.psi_old{k}(ind_m);
                    c_adi    {k} = obj.psi    {k}(ind_m);
                end
                
                % Change from adiabatic to super-adi rep
                c_old = ham.change_rep( c_adi_old , obj.u_sa_old(:,:,ind_m) );
                c     = ham.change_rep( c_adi     , obj.u_sa    (:,:,ind_m) );
                
                rho_old_m = abs(c_old{m}).^2;
                rho_old_n = abs(c_old{n}).^2;
                rho_n     = abs(c{n}    ).^2;
                
                p = ( rho_n - rho_old_n ) ./ rho_old_m;
                
                probable = max( p , probable );
                
            end
            
        end
        
        % Preprocessing: before hopping
        function prep_hop ( obj,first_call )
            global time 

            prep_hop@sht.mssh( obj,first_call )
            
            % If eq_motion = super_adi then u_sa has already been computed
            if ~isa(time.eq_motion,'tmp.traj.eq_motion.super_adi')
                obj.u_sa = ham.pot_eig_sa(obj.mom, obj.pot_mat, obj.frc_mat, ...
                                    obj.D_new, obj.U_new);
            end
            
        end
        
    end
end

