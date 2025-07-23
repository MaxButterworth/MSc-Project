% *************************************************************************
%
% Linear Vibronic Coupling model in 2 dimensions
%
% See: On the inclusion of the diagonal Born-Oppenheimer correction in surface hopping methods
%  by: Rami Gherib; Liyuan Ye; Ilya G. Ryabinkin; Artur F. Izmaylov
%      J. Chem. Phys. 144, 154103 (2016)
%      https://doi.org/10.1063/1.4945817
%                                                  
% *************************************************************************

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2024 Burkhard Schmidt's group
%
% see the README file for license details.

classdef lvc_2D < pot.generic & handle
    
    properties (Access = public)
        
        W1
        W2
        A
        C
        D
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = lvc_2D
            obj.empty = false;   
            obj.W1 = 1;
            obj.W2 = 1;
            obj.A = 0;
            obj.C = 1;
            obj.D = 0;
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
                
                prt.disp ('Generic E x e conical intersection example')
                
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
            
            w1 = obj.W1;
            w2 = obj.W2;
            a  = obj.A;
            c  = obj.C;
            d  = obj.D;
            
            if obj.row==1 && obj.col==1
                V = 0.5 * ( w1.^2 .* ( x + a/2 ).^2 + w2.^2 .* y.^2 + d );
            elseif obj.row==2 && obj.col==2
                V = 0.5 * ( w1.^2 .* ( x - a/2 ).^2 + w2.^2 .* y.^2 - d );
            else
                V = c .* y;
            end
            
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            global hamilt space
            
            x = r{1};
            y = r{2};
            
            w1 = obj.W1;
            w2 = obj.W2;
            a  = obj.A;
            c  = obj.C;
            
            if obj.row==1 && obj.col==1
                F{1} = - ( w1.^2 .* (a + 2*x) ) / 2;
                F{2} = - w2.^2 .* y;
            elseif obj.row==2 && obj.col==2
                F{1} = ( w1.^2 .* (a - 2*x) ) / 2;
                F{2} = - w2.^2 .* y;
            else
                F{1} = 0;
                F{2} = - c;
            end
            
        end        
        
        % Evaluate derivative of the forces as negative gradients of potential
        function G = G(obj,r)
            global hamilt space
            
            x = r{1};
            
            w1 = obj.W1;
            w2 = obj.W2;
            
            G{2,1} = zeros(size(x));
            G{1,2} = zeros(size(x));
            
            if obj.row == obj.col
                G{1,1} = - w1.^2 * ones(size(x));
                G{2,2} = - w2.^2 * ones(size(x));
            else
                G{1,1} = zeros(size(x));
                G{2,2} = zeros(size(x));
            end
            
        end        
        
    end
end
