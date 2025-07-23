% *************************************************************************
%
% Spin-Boson Coupling
%
% See: Surface-hopping dynamics of a spin-boson system
%  by: Donal Mac Kernan; Giovanni Ciccotti; Raymond Kapral
%      J. Chem. Phys. 116, 2346–2353 (2002)
%      https://doi.org/10.1063/1.1433502
%                                                  
% *************************************************************************

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2024 Burkhard Schmidt's group
%
% see the README file for license details.

classdef spin_boson < pot.generic & handle
    
    properties (Access = public)
        
        Omega
        W_max
        Z
        
        w_0
        w_vec
        c_vec
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = spin_boson
            obj.empty = false;   
            obj.Omega = 1/3;
            obj.W_max = 3;
            obj.Z     = 0.1;
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            
            global hamilt space
            
            if hamilt.coupling.n_eqs ~= 2
                prt.error ('This potential is only for 2 (coupled) channels')
            end

            obj.w_0 = ( 1 - exp( - obj.W_max ) ) / space.n_dim;

            obj.w_vec = zeros(space.n_dim,1);
            obj.c_vec = zeros(space.n_dim,1);
            for j = 1:space.n_dim
                obj.w_vec(j) = - log( 1 - j * obj.w_0 );
                obj.c_vec(j) = sqrt( obj.Z * obj.w_0 ) * obj.w_vec(j);
            end
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            
            disp@pot.generic (obj)
            if obj.row==1 && obj.col==1
                
                prt.disp ('TBD')
                
            else
                prt.disp ('Same as above')
                prt.disp ('***************************************************************')
            end
            prt.disp (' ')
            
        end
        
        % Evaluate potential energy functions
        function V = V(obj,r)
            global space

            V = zeros(size(r{1}));
            
            Omg     = obj.Omega;
            w       = obj.w_vec;
            c       = obj.c_vec;
            
            if obj.row==1 && obj.col==1
                for j = 1:space.n_dim
                    V = V + 0.5 * w(j)^2 .* r{j}.^2 - c(j) .* r{j};
                end
            elseif obj.row==2 && obj.col==2
                for j = 1:space.n_dim
                    V = V + 0.5 * w(j)^2 .* r{j}.^2 + c(j) .* r{j};
                end
            else
                V = -Omg;
            end
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            global space

            F = cell(space.n_dim,1);

            w = obj.w_vec;
            c = obj.c_vec;
            
            if obj.row==1 && obj.col==1
                for j = 1:space.n_dim
                    F{j} = - w(j)^2 .* r{j} + c(j);
                end
            elseif obj.row==2 && obj.col==2
                for j = 1:space.n_dim
                    F{j} = - w(j)^2 .* r{j} - c(j);
                end
            else
                for j = 1:space.n_dim
                    F{j} = zeros(size(r{1}));
                end
            end
            
        end        
        
        % Evaluate derivatives of the forces
        function G = G(obj,r)
            global space

            G = cell(space.n_dim,space.n_dim);

            for i = 1:space.n_dim
                
                for j = i:space.n_dim
                    G{i,j} = zeros(size(r{1}));
                    G{j,i} = zeros(size(r{1}));
                end

                if obj.row==obj.col
                    G{i,i} = - obj.w_vec(i)^2;
                end

            end
            
        end        
        
    end
end
