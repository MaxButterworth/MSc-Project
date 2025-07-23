%--------------------------------------------------------------------------
%
% Compute diabatic force matrix at positions r
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function frc_dia_m = frc_dia_m(pos, m)
global hamilt space

% Calculate force in diabatic representation
if hamilt.pot{m,m}.empty % free particle
    frc_dia_m = cell(space.n_dim,1);
    for d=1:space.n_dim
        frc_dia_m{d} = zeros(size(pos{1}));
    end
else
    frc_dia_m = F ( hamilt.pot{m,m}, pos );
end
end
        
