%% Time-independent  Schrödinger equation: solve directly
cd TISE/Excited;
qm_setup('wave'); qm_init(3); qm_bound(); qm_cleanup();
cd ../..;
%% Time-independent  Schrödinger equation: propagate in imginary time
cd TISE/Cyclic;
qm_setup('wave'); qm_init(6); qm_propa('cheby_imag'); qm_cleanup();
cd ../..;
cd TISE/Linear;
qm_setup('wave'); qm_init(6); qm_propa('cheby_imag'); qm_cleanup();
cd ../..;
%% Time-dependent Schrödinger equation: propagate in real time
cd TDSE/Cyclic;
qm_setup('wave'); qm_init(6); qm_propa('strang'); qm_cleanup();
cd ../..;