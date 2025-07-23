%--------------------------------------------------------------------------
%
% Potential energy function:
% Gaussian bell-shaped function 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2020-20xy Burkhard Schmidt
%
% see the README file for license details.

classdef gauss < pot.generic & handle
    
    properties (Access = public)
        height      % Height of Gaussian
        pos_0       % Central position of Gaussian
        width       % Width of Gaussian
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = gauss
            obj.empty = false;
            obj.height = 1;
            obj.pos_0  = 0;
            obj.width  = 1;
        end
        
        % Initialize potential: Set/check parameters
        function init (obj)
            
            global space
            
            if ~isscalar(obj.height)
                prt.error ('Incompatible dimensionality for height')
            end
            
            if length(obj.pos_0)~=space.n_dim
                prt.error ('Incompatible dimensionality for pos_0')
            end
            
            if length(obj.width)~=space.n_dim
                prt.error ('Incompatible dimensionality for width')
            end
                        
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic(obj)
            prt.disp ( 'Gaussian bell-shaped potential energy function in N-dim' )
            prt.disp ( '***************************************************************' )
            prt.disp ( '                                                      ' )
            prt.disp ( '               N      [   ( Ri-R0i )^2 ]              ' )
            prt.disp ( ' V (R) = H * Prod exp [ - (--------)   ]              ' )
            prt.disp ( '              i=1     [   (  2*Wi  )   ]              ' )
            prt.disp ( '                                                      ' )
            prt.disp ( 'where the product extends over all spatial dimensions ' )
            prt.disp ( '***************************************************************' )
            prt.disp (   ' ')
            prt.disp ( [ 'Height of Gaussian         H : ' num2str(obj.height) ] )
            prt.disp ( [ 'Mean value position       R0 : ' num2str(obj.pos_0 ) ] )
            prt.disp ( [ 'Position uncertainty       W : ' num2str(obj.width ) ] )
            prt.disp (' ')
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)
            global space
            
            V = obj.height * ones ( size(r{1}) );
            
            % Tensor product of one-dimensional Gaussians
            for k = 1:space.n_dim
                V = V .* exp (  -((r{k}-obj.pos_0(k)) / (2*obj.width(k))).^2  );
            end
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            
            F = cell(size(r));
            
            pot = V(obj,r);     
            
            for d=1:length(r)
                F{d} = pot .* ( r{d}-obj.pos_0(d) ) / ( 2 * obj.width(d)^2 );
            end
        end
        
        % Evaluate derivative of the forces, i.e. negative Hessian matrix
        function G = G(obj,r)
            global space
            G = cell(space.n_dim,space.n_dim);
            
            pot   = V(obj,r); 
            
            for i = 1:space.n_dim
                G{i,i} = 1 / ( 2 * obj.width(i)^2 ) .* pot ...
                         + ( r{i}-obj.pos_0(i) ).^2 / ( 2 * obj.width(i)^2 )^2 .* pot;
                for j = (i+1):space.n_dim
                    G{i,j} = pot .* ( r{i}-obj.pos_0(i) ) / ( 2 * obj.width(i)^2 ) ...
                                 .* ( r{j}-obj.pos_0(j) ) / ( 2 * obj.width(j)^2 );
                    G{j,i} = G{i,j};
                end
            end
        end      

    end
end


