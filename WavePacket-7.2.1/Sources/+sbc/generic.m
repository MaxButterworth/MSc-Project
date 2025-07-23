%--------------------------------------------------------------------------
%
% Generic properties of all system-bath coupling class definitions
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.


classdef generic < handle
    
    properties (Access = public) 
        row         % row index (in diabatic potential matrix)
        col         % column index (in diabatic potential matrix)
        dvr         % Grid representation (in N dimensions)
        empty       % Is this system-bath coupling empty?
    end
    
    methods (Access = public)

        % Constructor: Set empty system-bath coupling as the default
        function obj = generic
            obj.empty = true;
        end

        % Initialize: dummy method
        function init(~)
        end
        
        % Display system-bath coupling, overloading default disp method
        function disp (obj)
            global hamilt
            prt.disp('***************************************************************')
            if obj.row == obj.col
                prt.disp(['System-bath coupling for channel ' hamilt.coupling.labels{obj.row} ':'])
            else
                prt.disp(['System-bath coupling for channels ' hamilt.coupling.labels{obj.row} ', ' hamilt.coupling.labels{obj.col} ':'])
            end
            if obj.empty
                prt.disp ('Not available')
                prt.disp ('***************************************************************')
                prt.disp (' ')
            end
        end

        % Grid representation of system-bath coupling
        function grid (obj)
            global space
            
            if ~obj.empty
                obj.dvr = chi ( obj, space.dvr );
            end
            
        end
        
    end
    
end

