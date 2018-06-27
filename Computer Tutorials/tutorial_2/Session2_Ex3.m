%program which does Computer tutorial 2, Exercise 3
%Normal linear regression model with indep. Normal Gamma prior
%Gibbs sampling for independent Normal-Gamma prior
%load in the data set. Here use house price data from hprice.txt

load hprice.txt;
n=size(hprice,1);
y=hprice(:,1);
x=hprice(:,2:5);
x=[ones(n,1) x];
k=size(x,2);

%Hyperparameters for independent Normal-Gamma prior
v0=5;
b0=0*ones(k,1);
b0(2,1)=10;
b0(3,1)=5000;
b0(4,1)=10000;
b0(5,1)=10000;
s02=1/4.0e-8;
capv0=(10000^2)*eye(k);
capv0(2,2)=25;
capv0(3,3)=2500^2;
capv0(4,4)=5000^2;
capv0(5,5)=5000^2;
capv0inv=inv(capv0);

%Ordinary least squares quantities
bols = inv(x'*x)*x'*y;
s2 = (y-x*bols)'*(y-x*bols)/(n-k);
v=n-k;

%Calculate a few quantities outside the loop for later use
xsquare=x'*x;
v1=v0+n;
v0s02=v0*s02;


%Now start Gibbs loop
%beta conditional on h is Normal
%h conditional on beta is Normal

%store all draws in the following matrices
%initialize them here
b_=[];
h_=[];


%Specify the number of replications
%number of burnin replications
s0=1000;
%number of retained replications
s1=10000;
s=s0+s1;

%choose a starting value for h
hdraw=1/s2;  % initial guess

for i = 1:s
    %draw from beta conditional on h
    capv1inv = capv0inv+ hdraw*xsquare;
    capv1=inv(capv1inv);
    b1 = capv1*(capv0inv*b0 + hdraw*xsquare*bols);
    bdraw=b1 + norm_rnd(capv1); % norm_rnd from proffs website
     
    %draw from h conditional on beta
    s12 = ((y-x*bdraw)'*(y-x*bdraw)+v0s02)/v1;
    hdraw=gamrnd(.5*v1,2/(v1*s12));
   
    if i>s0
        %after discarding burnin, store all draws
        b_ = [b_ bdraw];
        h_ = [h_ hdraw];

    end
end


alldraws = [b_'];
%The function momentg is taken from LeSage's toolbox
%it inputs all Gibbs draws and produces posterior
%mean, standard deviation, nse and rne
%it calculates what the book calls S(0) in various ways
%see momentg.m for more details
result = momentg(alldraws);
means=[result.pmean]';
stdevs=[result.pstd]';
nse=[result.nse]';


%Print out whatever you want
'Hyperparameters for independent Normal-Gamma prior'
b0
capv0
v0
s02

'Posterior results based on Informative Prior'
'number of burnin replications'
s0
'number of included replications'
s1

'posterior mean, standard deviation '
'beta '
[means stdevs]


