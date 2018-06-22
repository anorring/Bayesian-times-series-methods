function y = Draw_Normal(mean,SIGMA)
% Generate random numbers from Normal distribution
%------------------------------------------------------
% Make sure you know what you are doing before you use this function
% anywhere else. This function does not check for wrong input (for example
% non-square covariance matrix SIGMA).
%------------------------------------------------------

p = size(SIGMA,1);
y = mean + chol(SIGMA)'*rand(p,1);