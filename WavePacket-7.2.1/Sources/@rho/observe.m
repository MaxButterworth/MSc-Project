%--------------------------------------------------------------------------
%
% For a given density matrix x and for a given system
% matrices A,B,N,C, this function calculates the
% all the linear(!) observables y
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2024 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function observe ( state, step )

% Note: x(t) is shifted with respect to equilibrium x_e
A = state.A;
C = state.C;
Q = state.Q;
x = state.x + state.x_equilib;
y = zeros (1,length(state.C));

% Linear observable: <c|x> may be complex
for len = 1:length(C) 
        if ~Q{len}
            y(len) = real(C{len}*x);
        else
            y(len) = abs(C{len}*x)^2;
        end
end
state.y (step,:) = y;

% Norm of rho "vector"
state.norm = norm (x); 

% Energy of state (without control field)
state.energy = dot( x, A * x);