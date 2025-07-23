%--------------------------------------------------------------------------
% Trigonometric (cosine-shaped) model for dipole moment
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%
% see the README file for license details.


classdef cosine < dip.generic & handle
    
    properties (Access = public)
        
        pre         % prefactor
        exp         % exponent
        mul         % multiplicity
        phi         % phase shift
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = cosine
            obj.empty = false;
            obj.pre = 1;
            obj.exp = 1;
            obj.mul = 1;
            obj.phi = 0;
        end
        
        % Initialize dipole moment: Set/check parameters
        function init (~)
        end
        
        % Display dipole moment, overloading default disp method
        function disp (obj)
            disp@dip.generic(obj)
            prt.disp (' mu(Theta) = f * cos^n  (m*Theta+phi)               ')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp (['Prefactor     f : ' num2str(obj.pre)])
            prt.disp (['Exponent      n : ' num2str(obj.exp)])
            prt.disp (['Multiplicity  m : ' num2str(obj.mul)])
            prt.disp (['Phase shift phi : ' num2str(obj.phi)])
            prt.disp (' ')
        end
        
        % Evaluate dipole moment
        function mu = mu(obj,r)
            global space
            if isa (space.dof{1}, 'dof.legendre')
                if obj.phi ~= 0 || obj.mul ~=1
                    prt.error('Code missing')
                end
                mu = obj.pre * r{1} .^obj.exp;
            else
                mu = obj.pre * cos(obj.mul*r{1}+obj.phi).^obj.exp;
            end
            
        end
    end
end
