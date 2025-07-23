%--------------------------------------------------------------------------
%
% For a given state vector x and for a given system
% matrices A,B,N,D, this function calculates the
% all the quadratic(!) observables y
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
D = state.D;
x = state.x + state.x_equilib;
y = zeros (1,length(D));

% Quadratic observable: real if D is Hermitian
for len = 1:length(D) 
    y(len) = dot( x, D{len} * x);
end
state.y (step,:) = y;

% Norm of ket vector
state.norm = norm (x); 

% Energy of state (without control field)
state.energy = 1i * dot( x, A * x);
