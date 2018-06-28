%This code, written by Joshua Chan, estimates the VAR-SV of Exercise 20.7
clear; clc;
p = 2;          % if p > 4, need to change Y0 and Y below
nsim = 20000;
burnin = 1000;

    % load data

[Yraw,Dates] = xlsread('US_macrodata1.csv');
Y0=Yraw(1:4,:); % save the first 4 obs as the initial conditions
Y = Yraw(5:end,:);

[T,n] = size(Y);
y = reshape(Y',T*n,1);
k = n*p+1;         % # of coefficients in each equation
m = n*(n-1)/2;     % # of free elements in L
    % prior
a0 = zeros(m,1); iVa = eye(m);
beta0 = zeros(n*k,1);
    % precision for coefficients = 1; for intercepts = 1/10
tmp = ones(k*n,1);  tmp(1:p*n+1:k*n) = 1/10;  
iVbeta = sparse(1:k*n,1:k*n,tmp);
nu_h = 3*ones(n,1); S_h = .1^2*ones(n,1);
b0 = ones(n,1); iB0 = eye(n)/10;

    % compute and define a few things 
tmpY = [Y0(end-p+1:end,:); Y];
X_tilde = zeros(T,n*p); 
for i=1:p
    X_tilde(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X_tilde = [ones(T,1) X_tilde];
X = SURform2(X_tilde,n); 
L_id = nonzeros(tril(reshape(1:n^2,n,n),-1));
L = eye(n); 
LiDL = zeros(T*n,n);

    % initialize for storage 
store_h = zeros(nsim,T,n); 
store_beta = zeros(nsim,k*n); 
store_sigh2 = zeros(nsim,n);

    % initialize the Markov chain 
beta = (X'*X)\(X'*y);
h0 = log(mean(reshape((y - X*beta).^2,n,T),2));
h = repmat(h0',T,1);
sigh2 = .1*ones(n,1);
a = zeros(m,1);
  
for isim = 1:nsim + burnin 
        % sample beta    
    L(L_id) = a;
    bigL = kron(speye(T),L);
    LiDL = bigL'*sparse(1:T*n,1:T*n,reshape(1./exp(h)',1,T*n))*bigL;
    XLiDL = X'*LiDL;
    Kbeta = iVbeta + XLiDL*X;
    beta_hat = Kbeta\(iVbeta*beta0 + XLiDL*y);
    beta =  beta_hat + chol(Kbeta,'lower')'\randn(n*k,1);
    
        % sample h  
    U = reshape(y - X*beta,n,T)';
    s = (L*U')';     
    ystar = log(s.^2 + .0001);
    for i=1:n        
        h(:,i) = SVRW(ystar(:,i),h(:,i),h0(i),sigh2(i));    
    end
    
        % sample a
    E = zeros(T*n,m);
    count_E = 0;
    for ii=1:n-1
        E(ii+1:n:end,count_E+1:count_E+ii) = -U(:,1:ii);        
        count_E = count_E+ii;
    end
    iD = sparse(1:T*n,1:T*n,reshape(1./exp(h)',1,T*n));
    Ka = iVa + E'*iD*E;
    a_hat = Ka\(iVa * a0 + E'*iD*reshape(U',T*n,1));   
    a = a_hat + chol(Ka,'lower')'\randn(m,1);    
    
        % sample sigh2 
    e2 = (h - [h0'; h(1:end-1,:)]).^2;             
    sigh2 = 1./gamrnd(nu_h + T/2, 1./(S_h + sum(e2)'/2));    
    
        % sample h0
    Kh0 = iB0 + sparse(1:n,1:n,1./sigh2);
    h0_hat = Kh0\(iB0*b0 + h(1,:)'./sigh2);
    h0 = h0_hat + chol(Kh0,'lower')'\randn(n,1);
    
    if (mod(isim, 5000) == 0)
        disp([num2str(isim) ' loops... ']);
    end    
      
        % store the parameters
    if isim > burnin
        isave = isim - burnin;
        store_beta(isave,:) = beta';
        store_h(isave,:,:) = h;
        store_sigh2(isave,:) = sigh2';       
    end
end
beta_hat = mean(store_beta)';
h_hat = squeeze(mean(exp(store_h/2)));
tt = linspace(1960,2013.75,T)';
figure;
for i=1:n
    subplot(1,n,i); plot(tt,h_hat(:,i));
    box off; xlim([1960 2014]);
end
set(gcf,'Position',[100 100 800 300]);