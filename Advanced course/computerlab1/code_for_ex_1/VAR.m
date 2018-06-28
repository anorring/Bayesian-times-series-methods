clear; clc;
p = 2;          % if p > 4, need to change Y0 and Y below
nsim = 2000;
burnin = 100;

    % load data

[Yraw,Dates] = xlsread('US_macrodata.csv');
Y0=Yraw(1:4,:); % save the first 4 obs as the initial conditions
Y = Yraw(5:end,:);

[T,n] = size(Y);
y = reshape(Y',T*n,1);
k = n*p+1;         % # of coefficients in each equation
n_hz = 20;         % # of horizons for IRs 
    % prior
nu0 = n+3; S0 = eye(n);
beta0 = zeros(n*k,1);
    % precision for coefficients = 1; for intercepts = 1/10
tmp = ones(k*n,1);  tmp(1:p*n+1:k*n) = 1/10;  
iVbeta = sparse(1:k*n,1:k*n,tmp);
    % compute and define a few things 
tmpY = [Y0(end-p+1:end,:); Y];
X_tilde = zeros(T,n*p); 
for i=1:p
    X_tilde(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X_tilde = [ones(T,1) X_tilde];
X = SURform2(X_tilde,n); 

    % initialize for storage 
store_Sig = zeros(nsim,n,n); 
store_beta = zeros(nsim,k*n); 
store_yIR = zeros(n_hz,n); 

    %% initialize the chain 
beta = (X'*X)\(X'*y);
e = reshape(y - X*beta,n,T);
Sig = e*e'/T;    
iSig = Sig\speye(n);
 
for isim = 1:nsim + burnin  
        % sample beta
    XiSig = X'*kron(speye(T),iSig);
    XiSigX = XiSig*X;
    Kbeta = iVbeta + XiSigX;
    beta_hat = Kbeta\(iVbeta*beta0 + XiSig*y); 
    beta =  beta_hat + chol(Kbeta,'lower')'\randn(n*k,1);

        % sample Sig
    e = reshape(y - X*beta,n,T);    
    Sig = iwishrnd(S0 + e*e',nu0 + T);
    iSig = Sig\speye(n);
    
    if (mod(isim, 1000) == 0)
        disp([num2str(isim) ' draws... ']);
    end
    
        % store the parameters
    if isim > burnin
        isave = isim - burnin;
        store_beta(isave,:) = beta';
        store_Sig(isave,:,:) = Sig;
        
            % compute impulse-responses
        CSig = chol(Sig,'lower');
                % 100 basis pts rather than 1 std. dev.
        shock = [0; 0; 1]/CSig(n,n); 
        yIR = construct_IR(beta,Sig,n_hz,shock);
        store_yIR = store_yIR + yIR;        
    end
end

beta_hat = mean(store_beta)';
Sig_hat = squeeze(mean(store_Sig));
yIR_hat = store_yIR/nsim;

figure;
subplot(1,3,1);
plot(yIR_hat(:,1)); box off; xlim([1 n_hz]);
subplot(1,3,2);
plot(yIR_hat(:,2)); box off; xlim([1 n_hz]);
subplot(1,3,3);
plot(yIR_hat(:,3)); box off; xlim([1 n_hz]);
set(gcf,'Position',[100 100 800 300]);
