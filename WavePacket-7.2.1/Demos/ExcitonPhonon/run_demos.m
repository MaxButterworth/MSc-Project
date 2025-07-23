%% Time-independent Schrödinger equation: direct diagonalization
cd TISE/Cyclic/N=03;
qm_setup('wave'); qm_init(3); qm_bound; qm_cleanup();
cd ../../..;
%% Time-independent Schrödinger equation: propagate in imginary time
cd TISE/Cyclic/N=06;
qm_setup('wave'); qm_init(6); qm_propa('cheby_imag'); qm_cleanup();
cd ../../..;
