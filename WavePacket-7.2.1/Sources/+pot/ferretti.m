% *************************************************************************
%
% Model of a conical intersection in two dimensions
%
% See: Quantum mechanical and semiclassical dynamics at a conical intersection
%  by: A. Ferretti; G. Granucci; A. Lami; M. Persico; G. Villani 
%      J. Chem. Phys. 104, 5517–5527 (1996)
%      https://doi.org/10.1063/1.471791
%
% *************************************************************************

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2024 Burkhard Schmidt's group
%
% see the README file for license details.

classdef ferretti < pot.generic & handle
    
    properties (Access = public)
        
        X1
        X2
        X3
        Kx
        Ky
        D
        
        Gamma
        A
        B
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = ferretti
            obj.empty = false;   
            obj.X1      = 0;
            obj.X2      = 0;
            obj.X3      = 0;
            obj.Kx      = 1;
            obj.Ky      = 1;
            obj.A       = 1;
            obj.B       = 1;
            obj.D       = 0;
            obj.Gamma   = 1;
        end
        
        % Initialize potential: Set/check parameters
        function init (~)
            
            global hamilt space
            
            if space.n_dim ~= 2
                prt.error ('This potential is only for 2 dimension')
            end
            
            if hamilt.coupling.n_eqs ~= 2
                prt.error ('This potential is only for 2 (coupled) channels')
            end
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            
            disp@pot.generic (obj)
            if obj.row==1 && obj.col==1
                
                prt.disp ('Model of a conical intersection in two dimensions')
                
            else
                prt.disp ('Same as above')
                prt.disp ('***************************************************************')
            end
            prt.disp (' ')
            
        end
        
        % Evaluate potential energy functions
        function V = V(obj,r)
            x = r{1};
            y = r{2};
            
            x1 = obj.X1;
            x2 = obj.X2;
            x3 = obj.X3;
            kx = obj.Kx;
            ky = obj.Ky;
            a  = obj.A;
            b  = obj.B;
            d  = obj.D;
            g  = obj.Gamma;
            
            if obj.row==1 && obj.col==1
                V = 0.5 * kx .* ( x - x2 ).^2 + 0.5 * ky .* y.^2 + d;
            elseif obj.row==2 && obj.col==2
                V = 0.5 * kx .* ( x - x1 ).^2 + 0.5 * ky .* y.^2;
            else
                V = g .* y .* exp( - a .* (x-x3).^2 - b .* y.^2 );
            end
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            x = r{1};
            y = r{2};
            
            x1 = obj.X1;
            x2 = obj.X2;
            x3 = obj.X3;
            kx = obj.Kx;
            ky = obj.Ky;
            a  = obj.A;
            b  = obj.B;
            g  = obj.Gamma;
            
            if obj.row==1 && obj.col==1
                F{1} = - kx .* (x - x2);
                F{2} = - ky .* y;
            elseif obj.row==2 && obj.col==2
                F{1} = - kx .* (x - x1);
                F{2} = - ky .* y;
            else
                expo = exp( - a .* (x-x3).^2 - b .* y.^2 );
                V12  = g .* y .* expo;
                
                F{1} = 2 * a .* (x-x3) .* V12;
                F{2} = - g .* expo + 2 * b .* y .* V12;
            end
            
        end        
        
        % Evaluate derivatives of the forces
        function G = G(obj,r)
            x = r{1};
            y = r{2};
            
            x3 = obj.X3;
            kx = obj.Kx;
            ky = obj.Ky;
            a  = obj.A;
            b  = obj.B;
            g  = obj.Gamma;
            
            G{2,1} = zeros(size(x));
            G{1,2} = zeros(size(x));
            
            if obj.row == obj.col
                G{1,1} = - kx * ones(size(x));
                G{2,2} = - ky * ones(size(x));
            else
                expo = exp( - a .* (x-x3).^2 - b .* y.^2 );
                V12  = g .* y .* expo;
                
                G{1,1} = 2 .* a .* V12 - 4 * a^2 .* (x-x3).^2 .* V12;
                G{2,2} = 6 .* b .* V12 - 4 * b^2 * y.^2 .* V12;
                
                G{1,2} = 2 * a * g * (x-x3) .* expo - 4 * a * b * (x-x3) .* y .* V12;
                G{2,1} = G{1,2};
            end
            
        end        
        
    end
end
