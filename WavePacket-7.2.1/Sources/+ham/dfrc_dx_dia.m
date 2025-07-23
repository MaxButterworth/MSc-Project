%--------------------------------------------------------------------------
%
% Compute derivative of diabatic force matrix at positions r
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function dfrc_dx_mats = dfrc_dx_dia(r)
global hamilt space

dfrc_dx_mats = cell(space.n_dim,space.n_dim);
for d=1:space.n_dim
    for k=1:space.n_dim
        dfrc_dx_mats{d,k} = zeros (hamilt.coupling.n_eqs,hamilt.coupling.n_eqs,length(r{1}));
    end
end

for p = 1:hamilt.coupling.n_eqs % diagonal
    for q = p:hamilt.coupling.n_eqs % off-diagonal
        G_pq = G ( hamilt.pot{p,q}, r );
        for d=1:space.n_dim
            for k=1:space.n_dim
                dfrc_dx_mats{d,k}(p,q,:) = G_pq{d,k};
                dfrc_dx_mats{d,k}(q,p,:) = G_pq{d,k};
            end
        end
    end
end
end
        
