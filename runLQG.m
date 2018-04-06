% Compute optimal controller and estimator for generalized LQG
%
% u(t)    = -L(t) x(t)
% x(t+1)  = lqgParams.A x(t) + B (I + Sum(lqgParams.C(i) rnd_1)) u(t) + C0 rnd_n
% y(t)    = H x(t) + Sum(D(i) rnd_1) x(t) + D0 rnd_n
% xhat(t+1) = lqgParams.A xhat(t) + B u(t) + K(t) (y(t) - H xhat(t)) + E0 rnd_n
% x(1)    ~ mean X1, covariance S1
%
% cost(t) = u(t)' R u(t) + x(t)' Q(t) x(t)
%
% NSim    number of simulated trajectories (default 0)  (optional)
% Init    0 - open loop; 1 (default) - LQG; 2 - random  (optional)
% Niter   iterations; 0 (default) - until convergence   (optional)
%
% K       Filter gains
% L       Control gains
% Cost    Expected cost (per iteration)
% Xa      Expected trajectory
% XSim    Simulated trajectories
% CostSim Empirical cost
% Unoise  Sequence of motor commands corrupted by noise
% Ufree   sequence of motor commands error FREE



% function [K,L,Cost,Xa,XSim,Xhat,CostSim,Unoise,Ufree] = ...
%    kalman_lqg(lqgParams.A,lqgParams.B,lqgParams.C,lqgParams.C0,lqgParams.H,lqgParams.D,lqgParams.D0, lqgParams.E0, lqgParams.Q, lqgParams.R, X1,lqgParams.S1, lqgParams.NSim,lqgParams.Init,lqgParams.Niter)

function [K,L,Cost,Xa,XSim,Xhat,CostSim,Unoise,Ufree] = runLQG(lqgParams, X1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization

% numerical parameters
MaxIter = 500;
Eps = 10^-15;

% determine sizes
szX = size(lqgParams.A,1);
szU = size(lqgParams.B,2);
szY = size(lqgParams.H,1);
szC = size(lqgParams.C,3);
szC0 = size(lqgParams.C0,2);
szD = size(lqgParams.D,3);
szD0 = size(lqgParams.D0,2);
szE0 = size(lqgParams.E0,2);
N = size(lqgParams.Q,3);

% initialize missing optional parameters
% if nargin<13,
%    lqgParams.NSim = 0;
% end;
% if nargin<14,
%    lqgParams.Init = 1;
% end;
% if nargin<15,
%    lqgParams.Niter = 0;
% end;

% if C or D are scalar, replicate them into vectors
if size(lqgParams.C,1)==1 & szU>1,
   lqgParams.C = lqgParams.C*ones(szU,1);
end;
if length(lqgParams.D(:))==1,
    if lqgParams.D(1)==0,
        lqgParams.D = zeros(szY,szX);
    else
        lqgParams.D = lqgParams.D*ones(szX,1);
        if szX ~= szY,
            error('D can only be a scalar when szX = szY');
        end
    end
end;

% if C0,D0,E0 are scalar, set them to 0 matrices and adjust size
if length(lqgParams.C0(:))==1 & lqgParams.C0(1)==0,
    lqgParams.C0 = zeros(szX,1);
end
if length(lqgParams.D0(:))==1 & lqgParams.D0(1)==0,
    lqgParams.D0 = zeros(szY,1);
end
if length(lqgParams.E0(:))==1 & lqgParams.E0(1)==0,
    lqgParams.E0 = zeros(szX,1);
end


% initialize policy and filter
K = zeros(szX,szY,N-1);
L = zeros(szU,szX,N-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run iterative algorithm - until convergence or MaxIter

for iter = 1:MaxIter
   
   % initialize covariances
   SiE = lqgParams.S1;
   SiX = X1*X1';
   SiXE = zeros(szX,szX);
   
   % forward pass - recompute Kalman filter   
   for k = 1:N-1
      
      % compute Kalman gain
      temp = SiE + SiX + SiXE + SiXE';
      if size(lqgParams.D,2)==1,
          DSiD = diag(diag(temp).*lqgParams.D.^2);
      else
          DSiD = zeros(szY,szY);
          for i=1:szD
              DSiD = DSiD + lqgParams.D(:,:,i)*temp*lqgParams.D(:,:,i)';
          end;
      end;
      K(:,:,k) = lqgParams.A*SiE*lqgParams.H'*pinv(lqgParams.H*SiE*lqgParams.H'+lqgParams.D0*lqgParams.D0'+DSiD);
      
      % compute new SiE
      newE = lqgParams.E0*lqgParams.E0' + lqgParams.C0*lqgParams.C0' + (lqgParams.A-K(:,:,k)*lqgParams.H)*SiE*lqgParams.A';
      LSiL = L(:,:,k)*SiX*L(:,:,k)';
      if size(lqgParams.C,2)==1,
         newE = newE + lqgParams.B*diag(diag(LSiL).*lqgParams.C.^2)*lqgParams.B';
      else
         for i=1:szC
            newE = newE + lqgParams.B*lqgParams.C(:,:,i)*LSiL*lqgParams.C(:,:,i)'*lqgParams.B';
         end;
      end;
      
      % update SiX, SiE, SiXE
      SiX = lqgParams.E0*lqgParams.E0' + K(:,:,k)*lqgParams.H*SiE*lqgParams.A' + (lqgParams.A-lqgParams.B*L(:,:,k))*SiX*(lqgParams.A-lqgParams.B*L(:,:,k))' + ...
          (lqgParams.A-lqgParams.B*L(:,:,k))*SiXE*lqgParams.H'*K(:,:,k)' + K(:,:,k)*lqgParams.H*SiXE'*(lqgParams.A-lqgParams.B*L(:,:,k))';
      SiE = newE;
      SiXE = (lqgParams.A-lqgParams.B*L(:,:,k))*SiXE*(lqgParams.A-K(:,:,k)*lqgParams.H)' - lqgParams.E0*lqgParams.E0';
   end;
   
   
   % first pass initialization
   if iter==1,
      if lqgParams.Init==0,         % open loop
         K = zeros(szX,szY,N-1); 
      elseif lqgParams.Init==2,     % random
         K = randn(szX,szY,N-1);
      end;
   end;
   
   
   % initialize optimal cost-to-go function
   Sx = lqgParams.Q(:,:,N);
   Se = zeros(szX,szX);
   Cost(iter) = 0;
   
   % backward pass - recompute control policy
   for k=N-1:-1:1
      
      % update Cost
      Cost(iter) = Cost(iter) + trace(Sx*lqgParams.C0*lqgParams.C0') + ...
         trace(Se*(K(:,:,k)*lqgParams.D0*lqgParams.D0'*K(:,:,k)' + lqgParams.E0*lqgParams.E0' + lqgParams.C0*lqgParams.C0'));
      
      % Controller
      temp = lqgParams.R + lqgParams.B'*Sx*lqgParams.B;
      BSxeB = lqgParams.B'*(Sx+Se)*lqgParams.B;
      if size(lqgParams.C,2)==1,
         temp = temp + diag(diag(BSxeB).*lqgParams.C.^2);
      else
         for i=1:size(lqgParams.C,3)
            temp = temp + lqgParams.C(:,:,i)'*BSxeB*lqgParams.C(:,:,i);
         end;
      end;
      L(:,:,k) = pinv(temp)*lqgParams.B'*Sx*lqgParams.A;
      
      % compute new Se
      newE = lqgParams.A'*Sx*lqgParams.B*L(:,:,k) + (lqgParams.A-K(:,:,k)*lqgParams.H)'*Se*(lqgParams.A-K(:,:,k)*lqgParams.H);
      
      % update Sx and Se
      Sx = lqgParams.Q(:,:,k) + lqgParams.A'*Sx*(lqgParams.A-lqgParams.B*L(:,:,k));
      KSeK = K(:,:,k)'*Se*K(:,:,k);
      if size(lqgParams.D,2)==1,
          Sx = Sx + diag(diag(KSeK).*lqgParams.D.^2);
      else
          for i=1:szD
              Sx = Sx + lqgParams.D(:,:,i)'*KSeK*lqgParams.D(:,:,i);
          end;
      end;     
      Se = newE;
   end;
   
   
   % adjust cost
   Cost(iter) = Cost(iter) + X1'*Sx*X1 + trace((Se+Sx)*lqgParams.S1);     
   
   % progress bar
   if ~rem(iter,10),
      fprintf('.');
   end;
   
   % check convergence of Cost
   if (lqgParams.Niter>0 & iter>=lqgParams.Niter) | ...
      (lqgParams.Niter==0 & iter>1 & abs(Cost(iter-1)-Cost(iter))<Eps) | ...
      (lqgParams.Niter==0 & iter>20 & sum(diff(dist(iter-10:iter))>0)>3),
      break;
   end;
end;

% % print result
% if Cost(iter-1)~=Cost(iter)
%    fprintf(' Log10DeltaCost = %.2f\n',log10(abs(Cost(iter-1)-Cost(iter))));
% else
%    fprintf(' DeltaCost = 0\n' );
% end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute average trajectory

Xa = zeros(szX,N);
Xa(:,1) = X1;

for k=1:N-1
   u = -L(:,:,k)*Xa(:,k);
   Xa(:,k+1) = lqgParams.A*Xa(:,k) + lqgParams.B*u;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate noisy trajectories

if lqgParams.NSim > 0,
   
   % square root of S1
   [u,s,v] = svd(lqgParams.S1);
   sqrtS = u*diag(sqrt(diag(s)))*v';
   
   % initialize
   XSim = zeros(szX,lqgParams.NSim,N);
   Xhat = zeros(szX,lqgParams.NSim,N);
   Xhat(:,:,1) = repmat(X1, [1 lqgParams.NSim]);
   XSim(:,:,1) = repmat(X1, [1 lqgParams.NSim]) + sqrtS*randn(szX,lqgParams.NSim);
   
   CostSim = 0;
   
   % loop over N
   for k=1:N-1
      
      % update control and cost
      U = -L(:,:,k)*Xhat(:,:,k);
      CostSim = CostSim + sum(sum(U.*(lqgParams.R*U)));
      CostU(k) = sum(sum(U.*(lqgParams.R*U)));
      
      CostSim = CostSim + sum(sum(XSim(:,:,k).*(lqgParams.Q(:,:,k)*XSim(:,:,k))));
      
      % compute noisy control
      Un = U;
      if size(lqgParams.C,2)==1,
         Un = Un + U.*randn(szU,lqgParams.NSim).*repmat(lqgParams.C,[1,lqgParams.NSim]);
      else
         for i=1:szC
            Un = Un + (lqgParams.C(:,:,i)*U).*repmat(randn(1,lqgParams.NSim),[szU 1]);
         end;
      end;
      
      % compute noisy observation
      y = lqgParams.H*XSim(:,:,k) + lqgParams.D0*randn(szD0,lqgParams.NSim);
      if size(lqgParams.D,2)==1,
         y = y + XSim(:,:,k).*randn(szY,lqgParams.NSim).*repmat(lqgParams.D,[1,lqgParams.NSim]);
      else
         for i=1:szD
            y = y + (lqgParams.D(:,:,i)*XSim(:,:,k)).*repmat(randn(1,lqgParams.NSim),[szY 1]);
         end;
      end;
      
      XSim(:,:,k+1) = lqgParams.A*XSim(:,:,k) + lqgParams.B*Un + lqgParams.C0*randn(szC0,lqgParams.NSim);
      Xhat(:,:,k+1) = lqgParams.A*Xhat(:,:,k) + lqgParams.B*U + K(:,:,k)*(y-lqgParams.H*Xhat(:,:,k)) + ...
          lqgParams.E0*randn(szE0,lqgParams.NSim);

      Unoise(:,:,k) = Un;
      Ufree(:,:,k)  = U;
   end;
   
   % final cost update
   CostSim = CostSim + sum(sum(XSim(:,:,N).*(lqgParams.Q(:,:,N)*XSim(:,:,N))));
   CostSim = CostSim / lqgParams.NSim;
   
else
   XSim = [];
   CostSim = [];
end;
