% Exercise 5
% Monte Carlo Integration. Assume that we have R draws from the posterior
% of a parameter theta. Lets generate this parameter from a Normal density
% with posterior mean 1 and variance 4.
theta_mean=1;
theta_var = 4;
number_of_draws = 100000; % the number of draws we obtain will affect the precision
                                                  % of the approximation using Monte Carlo integration
                        
theta_post= theta_mean + sqrt(theta_var)*randn(number_of_draws,1);

% Lets estimate the expected value of the function theta^2 using MC integration
theta_sq = (theta_post).^2;
MC_integration = mean(theta_sq);

% Here is another way of doing exactly as above, but using a "for loop"

theta_sq1=0;
for i=1:number_of_draws
    theta_draw=theta_mean+ sqrt(theta_var)*randn(1,1);
    theta_sq1 = theta_sq1 + theta_draw^2;
end
theta_sq1=theta_sq1/number_of_draws;


disp('Monte Carlo integration')
disp(MC_integration)
disp('Monte Carlo integration (alternative method)')
disp(theta_sq1)


% Pr(theta_sq>2)=1-Pr(theta_sq<=2)
Pr=sum(theta_sq>2)/number_of_draws;

disp('Pr(theta_sq>2)')
disp(Pr)
