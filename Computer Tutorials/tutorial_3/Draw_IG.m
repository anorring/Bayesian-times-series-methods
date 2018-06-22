function y = Draw_IG(mu,lambda)
% Generate a draw from the Inverse Gaussian (Wald) distribution
%----------------------------------------------------------------
% This method is from
% John R. Michael, William R. Schucany, and Roy W. Haas (1976) Generating 
% Random Variates Using Transformations with Multiple Roots. The American 
% Statistician 30: 88-90.
%
% Written by Dimitris Korobilis
% Universite Catholique de Louvain
% 23 October 2010

% First generate a chi_square(1) variate
v0 = randn.^2;

% Now find the two roots in Eq. (5) of the paper
x1 = mu + (.5*(mu.^2)*v0)/lambda - (.5*mu/lambda)*sqrt(4*mu*lambda*v0 + (mu.^2)*(v0.^2));
x2 = (mu.^2)/x1;

% Accept x1 with probability p1_v0, otherwise accept x2
p1_v0 = mu/(mu+x1);
if rand > p1_v0
    y = x1;
else
    y = x2;
end