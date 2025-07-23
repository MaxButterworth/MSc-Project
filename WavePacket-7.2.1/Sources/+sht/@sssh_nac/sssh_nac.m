%--------------------------------------------------------------------------
%
% "Single switch surface hopping" class definition
%
% Hopping probabilities are evaluated:
% at crossings of diabatic potentials OR
% at minima of gaps between adiabatic potentials
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

classdef sssh_nac < sht.sssh & handle
    
    properties (Access = public)
        nac
        nac_old       % Adiabatic energies (eigenvalues) of the previous step
        nac_oldold    % Adiabatic energies (eigenvalues) of the prev prev step
    end
    
    methods (Access = public)
        
        % Constructor: Set default values for class properties
        function obj = sssh_nac(lzv,n,seed)
            
            % Inherit from "generic" super class
            obj = obj@sht.sssh(lzv,n,seed);
                
            % Setting text strings
            obj.string0 = ['sssh_nac_' int2str(obj.lz_variant)];    
            obj.string4 = ['Single switch surface hopping: LZ variant NAC ' int2str(obj.lz_variant)];
        end
        
        function save_previous (obj) 
            
            % Inherit from initialization of superclass
            save_previous@sht.sssh ( obj);
            
            obj.nac_oldold  = obj.nac_old;
            obj.nac_old     = obj.nac;
        end
        
        % Preprocessing: before hopping
        function prep_hop ( obj,first_call )
            
            prep_hop@sht.sssh( obj,first_call )
            
            obj.nac = ham.nac_all (obj.pot_mat, obj.frc_mat, obj.U_new, obj.D_new);
            
        end
        
        function probable = prob_hop (obj,m,n,ind_m)
            global space time
            
            % Preallocate
            probable = zeros(size(ind_m));
            
            % Adiabatic energy gap
            if ~isempty(obj.nac_oldold)
                
                gap        = abs( obj.D_new(m,ind_m)   -obj.D_new(n,ind_m)    )';
                gap_old    = abs( obj.D_old(m,ind_m)   -obj.D_old(n,ind_m)    )';
                gap_oldold = abs( obj.D_oldold(m,ind_m)-obj.D_oldold(n,ind_m) )';
                
                nac_abs        = zeros(size(ind_m));
                nac_abs_old    = zeros(size(ind_m));
                nac_abs_oldold = zeros(size(ind_m));
                
                for d = 1:space.n_dim
                    nac_abs         = nac_abs        + abs(obj.nac       {m,n,d}(ind_m)).^2;
                    nac_abs_old     = nac_abs_old    + abs(obj.nac_old   {m,n,d}(ind_m)).^2;
                    nac_abs_oldold  = nac_abs_oldold + abs(obj.nac_oldold{m,n,d}(ind_m)).^2;
                end
                
                % Single switch citerion: Hopping can only occur at critical phase space points
                % Detect sign change of first derivative: from negative to positive
                ind_c = find ( nac_abs_oldold < nac_abs_old & nac_abs < nac_abs_old );
                
                if(~isempty(ind_c))
                    
                    pot_mats = obj.pot_mat(:,:,ind_m(ind_c));
                    
                    mom      = cell(space.n_dim,1);
                    frc_mats = cell(space.n_dim,1);
                    for d = 1:space.n_dim
                        mom{d}      = obj.mom{d}(ind_m(ind_c));
                        frc_mats{d} = obj.frc_mat{d}(:,:,ind_m(ind_c));
                    end
                    
                    % https://en.wikipedia.org/wiki/Finite_difference_coefficient
                    gap_d2  = (gap(ind_c) - 2* gap_old(ind_c) + gap_oldold(ind_c)) / time.steps.s_delta^2;
                    
                    probable (ind_c) = lz_formula(obj, mom,  pot_mats, frc_mats, ...
                                            obj.U_new(:,:,ind_m(ind_c)), obj.D_new(:,ind_m(ind_c)), ...
                                            gap(ind_c), gap_d2, m, n);
                end
            end
        end
    end
        
end
    
