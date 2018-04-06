function dnfParams=initDNFParams()

dnfParams.fieldSize = 181; % must be odd

dnfParams.tau_u = 20; % time constant of dynamic field
dnfParams.tau_mu_build = 500; dnfParams.tau_mu_decay = 2000; % time constants of memory traces
dnfParams.beta_u = 1; % steepness parameter of sigmoid function
dnfParams.h_u = -5; % resting level

% interation paramters
% c: strength of interaction kernel
% sigma: width of interaction kernel
% g: strength of global inhibition
dnfParams.c_exc = 10; dnfParams.sigma_exc = 5;
dnfParams.c_inh = 25; dnfParams.sigma_inh = 40; dnfParams.g_inh = 0.8;
dnfParams.c_mu = 0; dnfParams.sigma_mu = 5;

dnfParams.q_u = .25;    % noise levels
dnfParams.sigma_q = 5; % width of the noise kernel

dnfParams.halfField = floor(dnfParams.fieldSize/2);
