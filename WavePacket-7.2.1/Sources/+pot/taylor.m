%--------------------------------------------------------------------------
% Represent the potential energy function as a Taylor series
% Includes free particle (V=0) as a special case
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2016-2024 Burkhard Schmidt
%
% see the README file for license details.

classdef taylor < pot.generic & handle
    
    properties (Access = public)
        
        hshift      % Horizontal shift (row vector)
        vshift      % Vertical shift (scalar)
        coeffs      % Coefficients, i.e. derivatives
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = taylor
            obj.empty   = false;
            obj.hshift  = [];
            obj.vshift  = 0;
            obj.coeffs  = [];
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            global space
            if isempty(obj.hshift)
                obj.hshift = zeros(1,space.n_dim);
            end
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            prt.disp ('Taylor series (diag. in N dimensions)')
            prt.disp ('***************************************************************')
            prt.disp (' ')
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)  
            V = math.taylor (...
                r, ...
                obj.hshift, ...
                obj.vshift, ...
                obj.coeffs, 0 );
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            F = math.taylor_d (...
                r, ...
                obj.hshift, ...
                obj.coeffs );
        end
        
        % Evaluate derivative of the forces, i.e. negative Hessian matrix
        function G = G(obj,r)
            global space
            G = cell(space.n_dim,space.n_dim);
            [nrow,ncol]=size(obj.coeffs);
            
            for i = 1:space.n_dim
                G{i,i} = zeros(size(r{1}));
                for j = (i+1):space.n_dim
                    G{i,j} = zeros(size(r{1}));
                    G{j,i} = zeros(size(r{1}));
                end
            end
            
            for j=2:nrow
                fj = factorial(j);
                
                % Summing up contributions from each component of position vector
                for k = 1:ncol
                    G{k,k} = G{k,k} - j * (j-2) * (r{k}-obj.hshift(k)).^(j-2) * obj.coeffs(j,k) / fj;
                end
            end
            
        end
        
    end
end
