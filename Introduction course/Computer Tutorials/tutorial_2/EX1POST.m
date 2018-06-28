%For a given data set and values for prior hyperparameters
%Calculate OLS and posterior quantities, marginal likelihood
%and predictive (Predictive sets x-star=.5
%
% See lecture slides 27 -->

%Ordinary least squares quantities
bols = inv(x'*x)*x'*y;
s2 = (y-x*bols)'*(y-x*bols)/(n-k);
bolscov = s2*inv(x'*x);
bolssd=zeros(k,1);
for i = 1:k
bolssd(i,1)=sqrt(bolscov(i,i));
end
v=n-k;

%Posterior hyperparameters for Normal-Gamma
xsquare=x'*x;
v1=v0+n;
capv1inv = capv0inv+ xsquare;
capv1=inv(capv1inv);
b1 = capv1*(capv0inv*b0 + xsquare*bols);
if det(capv0inv)>0
    v1s12 = v0*s02 + v*s2 + (bols-b0)'*inv(capv0 + inv(xsquare))*(bols-b0);
else
    v1s12 = v0*s02 + v*s2;
end
s12 = v1s12/v1;

bcov = capv1*v1s12/(v1-2);
bsd=zeros(k,1);
for i = 1:k
bsd(i,1)=sqrt(bcov(i,i));
end

