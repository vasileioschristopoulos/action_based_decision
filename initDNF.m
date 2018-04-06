function dnf=initDNF(dnfParams, tStoreFields)
dnf.params=dnfParams;

% create row vectors for field activities
dnf.field_u = zeros(1, dnfParams.fieldSize);
dnf.memTrace_u = zeros(1, dnfParams.fieldSize);
dnf.output_u = zeros(1, dnfParams.fieldSize);

% create matrices to store field activities at different times
dnf.history_s = zeros(length(tStoreFields), dnfParams.fieldSize);
dnf.history_u = zeros(length(tStoreFields), dnfParams.fieldSize);
dnf.history_mu = zeros(length(tStoreFields), dnfParams.fieldSize);
dnf.history_output = zeros(length(tStoreFields), dnfParams.fieldSize);

% index of the current position in the history matrices
dnf.iHistory = 1;

% set up the interaction kernels
dnf.kSize_uu = min(round(3 * max(dnfParams.sigma_exc, dnfParams.sigma_inh)), floor((dnfParams.fieldSize-1)/2));
dnf.kernel_uu = dnfParams.c_exc * gaussNorm(-dnfParams.halfField:dnfParams.halfField, 0, dnfParams.sigma_exc) ...
    - dnfParams.c_inh * gaussNorm(-dnfParams.halfField:dnfParams.halfField, 0, dnfParams.sigma_inh) - dnfParams.g_inh;
dnf.kernel_mu = dnfParams.c_mu * gaussNorm(-dnfParams.halfField:dnfParams.halfField, 0, dnfParams.sigma_mu);

% set up the kernel for correlated noise (if required)
if dnfParams.sigma_q > 0
    dnf.kernel_q = gaussNorm(-dnfParams.halfField:dnfParams.halfField, 0, dnfParams.sigma_q);
end
