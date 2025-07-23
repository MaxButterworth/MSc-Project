% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Note that the 'cputime' command returns the 
%            total CPU time used by MATLAB.
%
% Note that the 'tic,toc' commands are a stopwatch timer
%           based on the wall-clock time in MATLAB.
%
% Copyright (C) 2004-2024 Burkhard Schmidt's group
%
% see the README file for license details.

function clock
global info

% CPU time elapsed (tt = total time)
info.tt = cputime-info.start_time;

% Date and time (nicely formatted)
dd = date;
c = clock;
hh = num2str(c(4),'%2.2i');
mm = num2str(c(5),'%2.2i');
ss = num2str(c(6),'%5.2f');

% Output (wt = wall time)
prt.disp ('***************************************************************');
if info.tt<0.1
    prt.disp ( 'Starting the clock ...' )
    tic
    info.wt = 0;
else
    prt.disp (['Elapsed seconds : ' num2str(info.tt,'%11.5e') ' seconds'] );
    info.wt = toc;
end
prt.disp ('***************************************************************');
prt.disp ( ' ' )
prt.disp (['Date            : ', dd] )
prt.disp (['Time            : ', hh, ':', mm, ':', ss] );
prt.disp ( ' ' )
