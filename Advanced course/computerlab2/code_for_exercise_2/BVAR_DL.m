%% Bayesian VAR with Dirichlet-Lapalce priors

clear; clc;
nloop = 25000;
burnin = 15000;

load USdata.csv % Load data
n = 3; % No. of endogenous variables
p = 1;  % No. of lags
y0 = USdata(1,:);
y = USdata(2:end,:);
T= size(y,1);
longy = reshape(y',T*n,1);
k = n+p*n^2;

%% Construct X
X = [ones(T,1) USdata(1:end-1,:)];
bigX = SURform2(X,n);

%% Priors
% beta
smalla = .5;

% Sigma
nu0 = n+3; 
S0 = speye(n);
newnu = T + nu0;

% Storage Matrix
store_beta = zeros(nloop-burnin,k);
store_invSi1 = zeros(nloop-burnin,n,n);
store_psi = zeros(nloop-burnin,k);
store_tau = zeros(nloop-burnin,1);
store_vartheta = zeros(nloop-burnin,k);
%% initialize the chain
Sig = cov(y);   
invSig = Sig\speye(n);
psi = 10*ones(k,1);
tau = 0.1;
vartheta = 10*ones(k,1);

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC.... ');
disp(' ' );
start_time = clock;

for loop = 1:nloop
%% Sample VAR
% sample beta
invVbeta = sparse(diag(1./(psi.*(vartheta.^2)*tau^2)));
XinvSig = bigX'*kron(speye(T),invSig);
XinvSigX = XinvSig*bigX;
invDbeta = invVbeta + XinvSigX;
cholbeta = chol(invDbeta,'lower')\speye(k);
betahat = (cholbeta'*cholbeta)*(XinvSig*longy);
beta = betahat + chol(invDbeta,'lower')'\randn(k,1); 

% sample Sig
err = reshape(longy - bigX*beta,n,T);
newS = S0 + err*err';
Sig = iwishrnd(newS,newnu);
invSig = Sig\speye(n);
% sample psi
nupsibeta = vartheta.*(tau./abs(beta));
invpsibeta = random('inversegaussian',nupsibeta,1);
 psi = 1./invpsibeta;
% sample tau
tau = gigrnd(k*(smalla-1),1,2*sum(abs(beta)./vartheta),1);
% sample vartheta
bigL = zeros(k,1);
for tv = 1:k
    bigL(tv,:) = gigrnd(smalla-1,1,2*abs(beta(tv,:)),1);
end
vartheta = bigL/sum(bigL);


%% Storage
if loop>burnin
    i = loop-burnin;
    store_beta(i,:) = beta';
    store_invSi1(i,:,:) = invSig;
    store_tau(i,:) = tau;
    store_psi(i,:) = psi';
    store_vartheta(i,:) = vartheta';

    
end
 if ( mod( loop, 2000 ) == 0 )
        disp(  [ num2str( loop ) ' draws... ' ] )
   end    
end
disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );

    
figure
subplot(4,3,1)
hist(store_beta(:,1))
title('\beta_1')
subplot(4,3,2)
hist(store_beta(:,2))
title('\beta_2')
subplot(4,3,3)
hist(store_beta(:,3))
title('\beta_3')
subplot(4,3,4)
hist(store_beta(:,4))
title('\beta_4')
subplot(4,3,5)
hist(store_beta(:,5))
title('\beta_5')
subplot(4,3,6)
hist(store_beta(:,6))
title('\beta_6')
subplot(4,3,7)
hist(store_beta(:,7))
title('\beta_7')
subplot(4,3,8)
hist(store_beta(:,8))
title('\beta_8')
subplot(4,3,9)
hist(store_beta(:,9))
title('\beta_9')
subplot(4,3,10)
hist(store_beta(:,10))
title('\beta_1_0')
subplot(4,3,11)
hist(store_beta(:,11))
title('\beta_1_1')
subplot(4,3,12)
hist(store_beta(:,12))
title('\beta_1_2')
    