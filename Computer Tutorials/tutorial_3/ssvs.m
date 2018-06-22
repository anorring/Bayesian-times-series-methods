%Stochstic Search Variable Selection program
%Uses cross country growth data from Fernandez, Ley and Steel (JAE, 2001)
%uses semi-default automatic choice for small large variances
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
xraw=rawdat(:,2:42);
bigk=size(xraw,2);
%bigk is the number of potential explanatory variables

%subtract mean from all regressors as in FLS
mxraw=mean(xraw);
sdxraw=std(xraw);
for i=1:bigk
    xraw(:,i)=(xraw(:,i) - mxraw(1,i))/sdxraw(1,i);
end

nobs=size(xraw,1);
xmat = [ones(nobs,1) xraw];

%Make the semi-automatic choices of "small" and "large" prior variances
%as recommended by George, Sun and Ni (2008, JOE)
%This uses OLS so will not work if bigk>nobs
xtxinv = inv(xmat'*xmat);
b_ols = xtxinv*xmat'*y;
s2_ols = (y-xmat*b_ols)'*(y-xmat*b_ols)/(nobs - bigk-1);
b_cov = s2_ols*xtxinv;
b_sd = sqrt(diag(b_cov));



%prior hyperparameters
V0= 10^2;
c1 = 0.1;
c2=10;
tau1 = c1*b_sd(2:bigk+1,1);
tau2 = c2*b_sd(2:bigk+1,1);
%prior hyperparameters for error precision
s_bar_2 = 0.01;
prior_dof = 0;
post_dof = prior_dof + nobs;

%starting values
gammas = ones(bigk,1);
sigsq = 0.0001;

%if SSVS=1 does SSVS in normal way
%if SSVS=0 imposes particular choice for gamma (e.g. based on initial run)
SSVS=1;
if SSVS==0
    load gammas.dat;
    for i=1:bigk
        if gammas(i,1)>0.5;
            gammas(i,1)=1;
        else
            gammas(i,1)=0;
        end
    end
end

iter = 1000;
burn = 100;
gammas_final = zeros(iter-burn,bigk);
betas_final = zeros(iter-burn,bigk+1);
sig_final = zeros(iter-burn,1);

for i = 1:iter;
    
    %-----------------
    %sample the betas
    %-----------------
    V_betas_temp = V0 ;
        for j=1:bigk;
            temp1=gammas(j)*tau2(j,1)^2 + (1-gammas(j))*tau1(j,1)^2;
    V_betas_temp = [V_betas_temp; temp1];
        end
    V_betas = diag(V_betas_temp);
    D_beta = inv(xmat'*xmat/sigsq + inv(V_betas));

    d_beta = xmat'*y/sigsq ;
    H_beta = chol(D_beta);
    betas = D_beta*d_beta + H_beta'*randn(bigk+1,1);
    
    %-----------------
    %sample the gammas
    %-----------------
    if SSVS==1
    for j = 1:bigk;
        numerator = normpdf(betas(j+1),0,tau2(j,1));
        denominator  = numerator + normpdf(betas(j+1),0,tau1(j,1));
        prob = numerator/ denominator;
        unif = rand(1,1);
        gammas(j,1) = .5*sign(prob-unif) + .5;
        
    end;
    end
    
    %---------------
    %sample the variance parameter
    %---------------
    %check the following, put in terms of Gamma
   
    post_sse = (y-xmat*betas)'*(y - xmat*betas)+ prior_dof*s_bar_2;
    post_shape = 2/post_sse;
    
    prec = gamrnd(0.5*post_dof,post_shape);
    sigsq=1/prec;
    
    if i > burn;
        betas_final(i-burn,:) = betas';
        sig_final(i-burn,1) = sigsq;
        gammas_final(i-burn,:) = gammas';
    end;
end;

bmean = mean(betas_final)';
b2mo = mean(betas_final.^2)';
bsd = sqrt(b2mo - bmean.^2);
sigmean=mean(sig_final);
gmean=mean(gammas_final)';
save gammas.dat gmean -ASCII;

'Mean of sigma^2'
sigmean
varcount=cumsum(ones(bigk,1));
'Posterior mean and standard deviation of beta follow by post mean gamma'
[varcount, bmean(2:bigk+1,1), bsd(2:bigk+1,1) gmean]


    