function dnf=resetDNF(dnf, tStoreFields)

% reset field activities to resting levels
dnf.field_u = zeros(1, dnf.params.fieldSize);
dnf.memTrace_u = zeros(1, dnf.params.fieldSize);
dnf.output_u = zeros(1, dnf.params.fieldSize);
dnf.iHistory = 1;

% create matrices to store field activities at different times
dnf.history_s = zeros(length(tStoreFields), dnf.params.fieldSize);
dnf.history_u = zeros(length(tStoreFields), dnf.params.fieldSize);
dnf.history_mu = zeros(length(tStoreFields), dnf.params.fieldSize);
dnf.history_output = zeros(length(tStoreFields), dnf.params.fieldSize);
