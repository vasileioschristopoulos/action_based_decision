function dnf=runDNF(dnf, stimulus, tStoreFields, t)

% calculation of field outputs
dnf.output_u = sigmoid(dnf.field_u, dnf.params.beta_u, 0);
        
% circular padding of outputs for convolution
dnf.output_u_padded = [dnf.output_u(dnf.params.halfField+2:dnf.params.fieldSize), dnf.output_u, ...
    dnf.output_u(:, 1:dnf.params.halfField)];
dnf.memTrace_u_padded = [dnf.memTrace_u(dnf.params.halfField+2:dnf.params.fieldSize), dnf.memTrace_u, ...
    dnf.memTrace_u(:, 1:dnf.params.halfField)];
        
% get endogenous input to fields by convolving outputs with interaction kernels
dnf.conv_uu = conv2(1, dnf.kernel_uu, dnf.output_u_padded, 'valid');
dnf.conv_mu = conv2(1, dnf.kernel_mu, dnf.memTrace_u_padded, 'valid');
        
% create field noise for this timestep
dnf.noise_u = dnf.params.q_u * randn(1, dnf.params.fieldSize);
if dnf.params.sigma_q > 0 % create spatially correlated noise by convolution
    dnf.noise_u_padded = [dnf.noise_u(dnf.params.halfField+2:dnf.params.fieldSize), dnf.noise_u, ...
        dnf.noise_u(:, 1:dnf.params.halfField)];
    dnf.noise_u = conv2(1, dnf.kernel_q, dnf.noise_u_padded, 'valid');
end
        
% update field activities
dnf.field_u = dnf.field_u + 1/dnf.params.tau_u * (-dnf.field_u + dnf.params.h_u + stimulus + ...
    dnf.conv_uu + dnf.conv_mu) + dnf.noise_u;
        
% update memory trace (only if there is activity in the field)
dnf.activeRegions_u = dnf.field_u > 0;
if any(dnf.activeRegions_u)
    dnf.memTrace_u = dnf.memTrace_u + 1/dnf.params.tau_mu_build * (-dnf.memTrace_u + dnf.output_u) .* ...
        dnf.activeRegions_u + 1/dnf.params.tau_mu_decay * (-dnf.memTrace_u) .* (1-dnf.activeRegions_u);
else
    dnf.memTrace_u = dnf.memTrace_u + 1/dnf.params.tau_mu_build * (-dnf.memTrace_u + dnf.output_u) + ...
        1/dnf.params.tau_mu_decay * (-dnf.memTrace_u);
end
        
% store field activities at the selected time steps
if any(tStoreFields == t)
    dnf.history_s(dnf.iHistory, :) = stimulus;
    dnf.history_u(dnf.iHistory, :) = dnf.field_u;
    dnf.history_mu(dnf.iHistory, :) = dnf.memTrace_u;
    dnf.history_output(dnf.iHistory, :) = dnf.output_u;
    dnf.iHistory = dnf.iHistory + 1;
end
