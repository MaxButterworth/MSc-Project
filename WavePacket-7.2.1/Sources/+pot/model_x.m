%--------------------------------------------------------------------------
%
% Subotnik
% pubs.acs.org/doi/10.1021/jp206557h
%
%--------------------------------------------------------------------------
%
% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2024 Burkhard Schmidt's group
%
% see the README file for license details.

classdef model_x < pot.generic & handle
    
    properties (Access = public)
        
        A           % Diabatic energy parameter
        B           % Diabatic range parameter
        C           % Coupling energy parameter
                
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = model_x
            obj.A = 0.03;
            obj.B = 1.6;
            obj.C = 0.005; 

            obj.empty = false;  
        end
        
        % Initialize potential: Set/check parameters
        function init (~)
            
            global hamilt space
            
            if space.n_dim ~= 1
                prt.error ('This potential is only for 1 dimension')
            end
            
            if hamilt.coupling.n_eqs ~= 3
                prt.error ('This potential is only for 3 states')
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
            x = r{1};

            m = min(obj.row,obj.col);
            n = max(obj.row,obj.col);

            a     = obj.A;
            b     = obj.B;
            c     = obj.C;

            if m==1 && n==1
                V =   a * ( tanh(b*x)     + tanh(b*(x+7)) );
            elseif m==2 && n==2
                V = - a * ( tanh(b*x)     + tanh(b*(x-7)) );
            elseif m==3 && n==3
                V = - a * ( tanh(b*(x+7)) - tanh(b*(x-7)) );
            elseif m==1 && n==2
                V = c * exp( -    x.^2 );
            elseif m==1 && n==3
                V = c * exp( -(x+7).^2 );
            elseif m==2 && n==3
                V = c * exp( -(x-7).^2 );
            end

        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            x = r{1};

            m = min(obj.row,obj.col);
            n = max(obj.row,obj.col);

            a     = obj.A;
            b     = obj.B;
            c     = obj.C;

            if m==1 && n==1
                F{1} = - a * b * ( 2 - tanh(b*x).^2     - tanh(b*(x+7)).^2 );
            elseif m==2 && n==2
                F{1} =   a * b * ( 2 - tanh(b*x).^2     - tanh(b*(x-7)).^2 );
            elseif m==3 && n==3
                F{1} =   a * b * (   - tanh(b*(x+7)).^2 + tanh(b*(x-7)).^2 );
            elseif m==1 && n==2
                F{1} = 2 * c *   x   .* exp( -    x.^2 );
            elseif m==1 && n==3
                F{1} = 2 * c * (x+7) .* exp( -(x+7).^2 );
            elseif m==2 && n==3
                F{1} = 2 * c * (x-7) .* exp( -(x-7).^2 );
            end

        end

    end
    
end


