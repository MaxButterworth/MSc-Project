%--------------------------------------------------------------------------
%
% Change a vector from rep 1 to rep 2 by using a transformation matrix U_12
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu Araujo
%
% see the README file for license details.

function c_rep_2 = change_rep( c_rep_1 , U_1_to_2 )
global hamilt
c_rep_2 = cell(size(c_rep_1));
for k = 1 : hamilt.coupling.n_eqs
    c_k = zeros(size(c_rep_1{k}));
    for l = 1 : hamilt.coupling.n_eqs
        c_k(:) = c_k + conj( squeeze( U_1_to_2(l,k,:) ) ) .* c_rep_1{l};
    end
    c_rep_2{k} = c_k;
end
end
