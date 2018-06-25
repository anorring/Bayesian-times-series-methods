% Exercise 3
%---------------------------------------------------------
% Code that calculates the posterior of an unknown 
% normal mean (variance is known). With asterisks ********
% are the pieces of the code that you should change at
% a first stage (of course, you may want to change the
% DGP of y as well & experiment with different parameters)
%---------------------------------------------------------
clear all; clc; %close all;
% Assume that we have generated data. Lets generate y from a N(3.5,4) density:
cons = 3.5; % cons is known, save it as a MATLAB variable for later use  
sigma = 2;  % sigma is known, save it as a MATLAB variable for later use   

% How many observations to generate???
obs = 1000;          %******** A: the more observation you have, the more precise the estimation and the prior will be

% Now generate the sample y from the normal density
y = cons + sigma*randn(obs,1);
       
% Step 1, set values for the prior hyperparameters.
% An 'uninformative' choice would be mu = 0 and tau = 100000.
% However experiment with more extreme values. What would happen if you
% impose mu to be much different than the true value (3.5) and with set a
% small prior variance? 
mu = -5;  % this is the prior mean             %********
tau = 1;  % this is the prior variance         %********
       
% Step 2. Now give an expression for the posterior of the parameter theta   
theta_mean = ((sigma^2 + tau^2)^(-1) )*((sigma^2)*mu + (tau^2)*mean(y));
theta_variance = ((sigma^2 + tau^2)^(-1) )*(sigma^2)*(tau^2);
   
% Knowing the posterior mean and variance are enough in order to characterize 
% the posterior of the parameters, since this is a Normal density. However we can easily    
% generate samples from the posterior and plot the full density
samp_no = 100000;  % define number of samples to obtain           %**** 
theta_posterior = theta_mean + sqrt(theta_variance)*randn(samp_no,1); 
   
% Now get an idea of how your posterior looks like
figure
hist(theta_posterior,30)
title('histogram of the samples from the parameter posterior')

% You can also get samples from your prior and plot those as well   
theta_prior = mu + tau*randn(samp_no,1);
figure
hist(theta_prior,30)
title('histogram of the samples from the parameter prior')


% You can also view how close your initial guess (the prior) is
% compared to the posterior
figure
dens1=normpdf([-20:0.001:20],mu,tau);
dens2=normpdf([-20:0.001:20],theta_mean,sqrt(theta_variance));
dens3=normpdf([-20:0.001:20],cons,sqrt(sigma));
plot([-20:0.001:20],dens1,'-',[-20:0.001:20],dens2,'--',[-20:0.001:20],dens3,'-.')
legend('prior density','posterior density', 'likelihood (data) density')
title('prior and posterior')

% Results for simulation exercise:
disp('True mean of theta is 3.5')
disp('    ')
disp(['Average of the posterior draws of theta is  ' num2str(mean(theta_posterior))])
disp('    ')
disp(['The MLE quantity is  ' num2str(mean(y))])