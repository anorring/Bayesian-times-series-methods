%LASSO
%Uses cross country growth data from Fernandez, Ley and Steel (JAE, 2001)
%Warning: the Gamma distribution is parameterized in different ways and
%Matlab's gamrnd, my textbook and the Lasso literature use different ones
%This accounts for the different ways things are drawn in this code
%see wikipedia's entry on Gamma distribution for more explanation
load growth.dat;
%The data set is arranged with data for each country taking up 6 lines
%The following makes it into an N by K matrix
n=72;
rawdat=zeros(n,42);
j=1;
for i=1:n
    rawdat(i,:)= [growth(j,:) growth(j+1,:) growth(j+2,:) ...
            growth(j+3,:) growth(j+4,:) growth(j+5,:)];
    j=j+6;
end

y=rawdat(:,1);
y=(y-mean(y));
xraw=rawdat(:,2:42);
bigk=size(xraw,2);
%bigk is the number of potential explanatory variables

%subtract mean from all regressors as in FLS
mxraw=mean(xraw);
sdxraw=std(xraw);
for i=1:bigk
    xraw(:,i)= (xraw(:,i) - mxraw(1,i))/sdxraw(1,i);
end

x = [ones(n,1) xraw];
T = size(x,1); % time series observations
p = size(x,2); % maximum nubmer of predictors

% ----------------Gibbs related preliminaries
nsave = 1000;
nburn = 100;
ntot = nsave + nburn;

beta_draws = zeros(nsave,p);
tau2_draws = zeros(nsave,p);
lambda2_draws = zeros(nsave,1);
sigma2_draws = zeros(nsave,1);
% ----------------Set priors
%this code is written with Gamma parameterized in terms of shape and inverse scale
% lambda2 ~ Gamma(r,d)
%These prior hyperparameter choices are not identical to those in the
%handout so your results with them may differ a bit
%You may want to experiment with these to investigate prior sensitivity
r = 1;
delta =1;
% Get OLS quanities from the full model (only if sample size is large enough)

if T>p
    beta_OLS = inv(x'*x)'*(x'*y);
    SSE_OLS = (y - x*beta_OLS)'*(y - x*beta_OLS);
    sigma2_OLS = SSE_OLS/(T-(p-1));
    b_cov = sigma2_OLS*inv(x'*x);
b_sd = sqrt(diag(b_cov));

end

% Initialize parameters
beta = 0*ones(p,1); 
tau2 = 4*ones(p,1);
V_L = diag(tau2);
lambda2 = 1;
sigma2 = 0.001;
if T>p
    beta = beta_OLS;
lambda2=(p*sqrt(sigma2_OLS)/(sum(abs(beta_OLS)))).^2;
sigma2=sigma2_OLS;
end


%==========================================================================
%====================| GIBBS ITERATIONS START HERE |=======================
disp('Now you are running the Bayesian Lasso')
for irep = 1:ntot
   if mod(irep,5000)==0
        disp(irep)
    end
    
    % 1. Update beta from Normal
    A = inv(x'*x + inv(V_L));
    post_mean_beta = A*x'*y;
    post_var_beta = sigma2*A;
    beta = Draw_Normal(post_mean_beta,post_var_beta);
    
    % 2. Update tau2_j from Inverse Gaussian
    for j = 1:p
        a1 = (lambda2*sigma2)./(beta(j,1)^2);
        a2 = lambda2;
        tau2_inverse = Draw_IG(sqrt(a1),a2);
 %note: often need to add a very small constant to avoid matrix singulariry
        tau2(j,1) = 1/tau2_inverse + 1e-15;
    end
    
    % Now that we have the new estimate of tau2, update V_L (the prior
    % covariance matrix of beta)
      V_L = diag(tau2);
    
    % 3. Update lambda2 from Gamma
    b1 = p + r;
    b2 = 0.5*sum(tau2) + delta;
   lambda2=gamrnd(b1,1/b2);
 %   lambda2 = Draw_Gamma(b1,b2);
       
    % 4. Update sigma2 from Inverse Gamma
    c1 = (T-1+p)/2;
    PSI = (y-x*beta)'*(y-x*beta);
    c2 = 0.5*PSI + 0.5*(beta'/V_L)*beta;  
    sigma2 = Draw_iGamma(c1,c2);
 
    % Save draws
    if irep > nburn
        beta_draws(irep-nburn,:) = beta;
        tau2_draws(irep-nburn,:) = tau2;
        lambda2_draws(irep-nburn,:) = lambda2;
        sigma2_draws(irep-nburn,:) = sigma2;
    end
    
end


% Now print some results
clc;
bmean=mean(beta_draws)';
b2mo = mean(beta_draws.^2)';
bsd = sqrt(b2mo - bmean.^2);
taumean=mean(sqrt(tau2_draws),1)';
varcount=cumsum(ones(bigk,1));
'Posterior mean and standard deviation of beta and mean of tau-s'
[varcount, bmean(2:bigk+1,1), bsd(2:bigk+1,1), taumean(2:bigk+1,1)]

if T>p
    'OLS estimates and standard errors'
    [varcount, beta_OLS(2:bigk+1,1), b_sd(2:bigk+1,1)]
    
end
disp('    ')
disp('    ')
disp('    ')

disp('POSTERIOR MEAN OF SIGMA2')
disp('    ')
disp(mean(sigma2_draws))
disp('POSTERIOR MEAN OF LAMBDA^2')
disp('    ')
disp(mean(lambda2_draws))


