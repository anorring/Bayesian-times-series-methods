%This code, written by Joshua Chan, estimates the TVP-VAR of Exercise 20.6
clear; clc;
p = 2;          % if p > 4, need to change Y0 and Y below
nsim = 20000;
burnin = 1000;
% load data

[Yraw,Dates] = xlsread('US_macrodata.csv');
Y0=Yraw(1:4,:); % save the first 4 obs as the initial conditions
Y = Yraw(5:end,:);
 
[T,n] = size(Y);
y = reshape(Y',T*n,1);
k = n*p+1;         % # of coefficients in each equation
t0 = [61 181]';    % 1975Q1, 2005Q1
n_hz = 20;         % # of horizons for IRs 

    % prior
nu0 = n+3; S0 = eye(n);
a0 = zeros(n*k,1);
    % precision for coefficients = 1; for intercepts = 1/10
tmp = ones(k*n,1);  tmp(1:p*n+1:k*n) = 1/10;  
iB0 = sparse(1:k*n,1:k*n,tmp);
nu0q = 3*ones(n*k,1); 
S0q = .01^2*ones(n*k,1); S0q(1:k:end) = .1^2;

    % compute and define a few things 
tmpY = [Y0(end-p+1:end,:); Y];
X_tilde = zeros(T,n*p); 
for i=1:p
    X_tilde(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X = SURform([ones(n*T,1) kron(X_tilde,ones(n,1))]);
H = speye(T*k*n,T*k*n) - sparse(n*k+1:T*k*n,1:(T-1)*n*k,...
    ones(1,(T-1)*n*k),T*n*k,T*n*k);

    % initialize for storage 
store_Sig = zeros(nsim,n,n); 
store_Q = zeros(nsim,k*n); 
store_beta = zeros(nsim,T*k*n); 
store_yIR_75 = zeros(n_hz,n); 
store_yIR_05 = zeros(n_hz,n);
store_diff = zeros(nsim,n_hz,n);

    % initialize the Markov chain 
Z = SURform2([ones(T,1) X_tilde],n);
beta0 = (Z'*Z)\(Z'*y);
e = reshape(y - Z*beta0,n,T);
Sig = e*e'/T;
iSig = Sig\speye(n);
Q = .1*ones(k*n,1);
 
for isim = 1:nsim + burnin  
        % sample beta
    HiQH = H'*sparse(1:T*n*k,1:T*n*k,repmat(1./Q,T,1))*H;
    XiSig = X'*kron(speye(T),iSig);
    Kbeta = HiQH + XiSig*X;
    beta_hat = Kbeta\(HiQH*kron(ones(T,1),beta0) + XiSig*y);    
    beta =  beta_hat + chol(Kbeta,'lower')'\randn(T*n*k,1);

        % sample Sig
    e = reshape(y - X*beta,n,T);    
    Sig = iwishrnd(S0 + e*e',nu0 + T);
    iSig = Sig\speye(n);
    
        % sample Q
    e = reshape(beta - [beta0;beta(1:end-n*k)],n*k,T);
    Q = 1./gamrnd(nu0q + T/2,1./(S0q + sum(e.^2,2)/2));    
        
        % sample beta0
    Kbeta0 = iB0 + sparse(1:n*k,1:n*k,1./Q);
    beta0_hat = Kbeta0\(iB0*a0 + sparse(1:n*k,1:n*k,1./Q)*beta(1:n*k)); 
    beta0 =  beta0_hat + chol(Kbeta0,'lower')'\randn(n*k,1);
    
    if (mod(isim, 1000) == 0)
        disp([num2str(isim) ' loops... ']);
    end
    
        % store the parameters
    if isim > burnin
        isave = isim - burnin;
        store_beta(isave,:) = beta';
        store_Sig(isave,:,:) = Sig;
        store_Q(isave,:) = Q';
        
           % compute impulse-responses
        CSig = chol(Sig,'lower');
                % 100 basis pts rather than 1 std. dev.
        shock = [0; 0; 1]/CSig(n,n); 
        tmp_beta = reshape(beta,n*k,T);
        yIR_75 = construct_IR(tmp_beta(:,t0(1)),Sig,n_hz,shock);
        yIR_05 = construct_IR(tmp_beta(:,t0(2)),Sig,n_hz,shock);
        
        store_yIR_75 = store_yIR_75 + yIR_75;
        store_yIR_05 = store_yIR_05 + yIR_05;
        store_diff(isave,:,:) = yIR_05 - yIR_75;
    end
end
beta_hat = mean(store_beta)';
yIR75_hat = store_yIR_75/nsim;
yIR05_hat = store_yIR_05/nsim;
diff_hat = squeeze(mean(store_diff));
diff_CI = quantile(store_diff,[.05 .95]);

figure;
subplot(1,3,1);
plot([yIR75_hat(:,1) yIR05_hat(:,1)]); box off; xlim([1 n_hz]);
subplot(1,3,2);
plot([yIR75_hat(:,2) yIR05_hat(:,2)]); box off; xlim([1 n_hz]);
subplot(1,3,3);
plot([yIR75_hat(:,3) yIR05_hat(:,3)]); box off; xlim([1 n_hz]);
set(gcf,'Position',[100 100 800 300]);

figure;
for i = 1:n
    subplot(1,n,i);
    hold on
        plot(diff_hat(:,i));
        plot(diff_CI(1,:,i),'--k');
        plot(diff_CI(2,:,i),'--k');
    hold off    
    box off; xlim([1 n_hz]);
end
set(gcf,'Position',[100 100 800 300]);

