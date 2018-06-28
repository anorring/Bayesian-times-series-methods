%program which does Computer tutorial 2, Exercise 1
%Normal linear regression model with natural conjugate prior (analytical)
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

%Call script which carries actually does posterior analysis
EX1POST;


%Print out whatever you want
'Hyperparameters for informative natural conjugate prior'
b0
capv0
v0
s02

'Posterior results based on Informative Prior'
b1
bsd




%Hyperparameters for noninformative prior
v0=0;
capv0inv=0*eye(k); % diagonal matrix, is this right? Or should the 


%Call script which carries actually does posterior analysis
EX1POST;

%Print out whatever you want
'Posterior results based on Noninformative Prior'
b1
bsd

