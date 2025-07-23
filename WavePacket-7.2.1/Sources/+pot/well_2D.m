% *************************************************************************
%
% Non-separable two-dimensional potential wells
%
% See: Phase-corrected surface hopping: Correcting the phase evolution of the electronic wavefunction
%  by: Neil Shenvi; Joseph E. Subotnik; Weitao Yang
%      J. Chem. Phys. 135, 024101 (2011)
%      https://﻿doi.org/10.1063/1.3603447
%                                                  
% *************************************************************************

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2024 Burkhard Schmidt's group
%
% see the README file for license details.

classdef well_2D < pot.generic & handle
    
    properties (Access = public)
        
        A
        B
        C
        D
        E_0

    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = well_2D
            obj.empty = false;   
            obj.A   = 0.15;
            obj.B   = 0.14;
            obj.C   = 0.015;
            obj.D   = 0.06;
            obj.E_0 = 0.05;
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            
            global hamilt space

            if space.n_dim ~= 2
                prt.error ('This potential is only for 2 dof')
            end
            
            if hamilt.coupling.n_eqs ~= 2
                prt.error ('This potential is only for 2 (coupled) channels')
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

            x=r{1};
            y=r{2};
            
            a     = obj.A;
            b     = obj.B;
            c     = obj.C;
            d     = obj.D;
            e_0   = obj.E_0;
            
            if obj.row==1 && obj.col==1
                V = - e_0 * ones(size(r{1}));
            elseif obj.row==2 && obj.col==2
                V = - a * exp( - b * ( 0.75*(x+y).^2 + 0.25*(x-y).^2 ) );
            else
                V =   c * exp( - d * ( 0.25*(x+y).^2 + 0.75*(x-y).^2 ) );
            end
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            global space

            F = cell(space.n_dim,1);
                
            x=r{1};
            y=r{2};
            
            a     = obj.A;
            b     = obj.B;
            c     = obj.C;
            d     = obj.D;
            
            if obj.row==1 && obj.col==1
                F{1} = zeros(size(x));
                F{2} = zeros(size(x));
            elseif obj.row==2 && obj.col==2
                V = - a * exp( - b * ( 0.75*(x+y).^2 + 0.25*(x-y).^2 ) );

                F{1} = 2 * b * (x + 0.5*y) .* V;
                F{2} = 2 * b * (y + 0.5*x) .* V;
            else
                V = c * exp( - d * ( 0.25*(x+y).^2 + 0.75*(x-y).^2 ) );

                F{1} = 2 * d * (x - 0.5*y) .* V;
                F{2} = 2 * d * (y - 0.5*x) .* V;
            end
            
        end        
        
        % Evaluate derivatives of the forces
        function G = G(obj,r)
            global space

            G = cell(space.n_dim,space.n_dim);

            prt.error ('Code missing')

%             for i = 1:space.n_dim
%                 
%                 for j = i:space.n_dim
%                     G{i,j} = zeros(size(r{1}));
%                     G{j,i} = zeros(size(r{1}));
%                 end
% 
%                 if obj.row==obj.col
%                     G{i,i} = - obj.w_vec(i)^2;
%                 end
% 
%             end
            
        end        
        
    end
end
