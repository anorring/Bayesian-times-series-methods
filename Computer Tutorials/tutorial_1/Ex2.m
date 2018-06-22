    %This is a program which artificially creates a data set and then does OLS estimation using it
    %First part of this program artificially simulates data
    n=100;
    alpha=1;
    beta=2;
    e = randn(n,1);
    x=rand(n,1);
    y=alpha + x*beta + e;
    %following line adds intercept to x. explain why
    x=[ones(n,1), x];
    %Following part of the program does OLS estimation
    %the OLS estimator
    bhat = inv(x'*x)*x'*y;
    disp 'The OLS estimator of beta is';
    disp(bhat);
    %the OLS residuals
    resids = y - x*bhat;
    %The OLS estimator of the error variance
    s2 = resids'*resids/(n-2);
    disp 'The OLS estimator of the error variance is';
    disp(s2);