%program which does Computer tutorial 2, Exercise 2
%Normal linear regression model with natural conjugate prior (Monte Carlo)
%load in the data set. Here use house price data from hprice.txt

load hprice.txt;
n=size(hprice,1);
y=hprice(:,1);
x=hprice(:,2:5);
x=[ones(n,1) x];
k=size(x,2);

%Hyperparameters for natural conjugate prior
v0=5;
b0=0*ones(k,1);
b0(2,1)=10;
b0(3,1)=5000;
b0(4,1)=10000;
b0(5,1)=10000;
s02=1/4.0e-8;
capv0=2.4*eye(k);
capv0(2,2)=6e-7;
capv0(3,3)=.15;
capv0(4,4)=.6;
capv0(5,5)=.6;
capv0inv=inv(capv0);
%Ordinary least squares quantities
bols = inv(x'*x)*x'*y;
s2 = (y-x*bols)'*(y-x*bols)/(n-k);
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

vscale = s12*capv1;
vchol=chol(vscale);
vchol=vchol';

%Now start Monte Carlo loop
%beta is t(b1,vscale,v1)
%For illustrative purpose we calculate only for beta(2)
 b2mean=zeros(k,1);
b2square=zeros(k,1);

%Specify the number of replications
s=100000;

%tdis_rnd takes random draws from the multivariate t
%with mean zero and scale, V=I
%Hence we have to transform draws to get t(b1,bscale,v1)
for i = 1:s
    %draw a t(v1) then transform to yield draw of beta
    
    
    bdraw=b1 + vchol*trnd(v1,k,1);
    b2mean=b2mean+bdraw;  
    b2square=b2square+bdraw.^2; 
end

b2mean=b2mean./s;
b2square=b2square./s;
b2var=b2square - b2mean.^2;
b2sd=sqrt(b2var);



%Print out whatever you want
'Hyperparameters for informative natural conjugate prior'
b0
capv0
v0
s02

'Posterior results based on Informative Prior'
b1
bsd
s
b2mean
b2sd
