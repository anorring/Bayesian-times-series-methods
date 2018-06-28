%program which does empirical illustration for the TVP-AR model

load dliprod.dat;
t=size(dliprod,1);
p=1;
dliprod=100*dliprod;
y=dliprod(p+1:t,1);
xar=ones(t-p,p+1);
for i=1:p
    xar(:,i+1) = dliprod(p+1-i:t-i,1);
end
t=t-p;
%k=total number of time varying parameters
k=p+1;

%specify hyperparameters of Gamma prior for h
v0=1;
s02=1;
%specify hyperparameters for Gamma prior for the vector lambda
vlam0=ones(k,1);
s02lam=1*ones(k,1);

%Specify some of the fixed matrices in the state space model
v1=v0+t;
vlam1=vlam0+t;

Gtt =zeros(1,k+1);
Gtt(1,1)=1;
Zt=xar;
Tt=eye(k);
%a temporary choice
Ht = zeros(k,k+1);
Ht(:,2:k+1)=eye(k);
Ft=Ht;


%store all draws in the following matrices
%initialize them here
h_=[];
lam_=[];

%Specify the number of replications
%number of burnin replications
s0=100;
%number of retained replications
s1=1000;
s=s0+s1;

%choose a starting value for h
hdraw=1;
adraw=zeros(t,k);
lamdraw=zeros(k,1);

for irep = 1:s
    irep
    %draw states conditional on other parameters
    sig2draw=1/hdraw;

    ndraw=djs(y,Zt,Gtt,Tt,Ht,Ft,sqrt(sig2draw));
    adraw(1,:)=ndraw(1,:);
    for i=2:t
      adraw(i,:)=adraw(i-1,:)+ ndraw(i,:);
    end
     
    %draw from h conditional on states
    sse=0;
    for i=1:t
        sse=sse+ (y(i,1) - xar(i,:)*adraw(i,:)')^2;
    end
    s12 = (sse+v0*s02)/v1;
    hdraw=gamm_rnd(1,1,.5*v1,.5*v1*s12);
    
    %draw from lambda conditional on rest
    sselam=(adraw(1,:)').^2;
    for i=2:t
      sselam=sselam + (adraw(i,:)' - adraw(i-1,:)').^2;
    end
    
    sselam=vlam0.*s02lam + hdraw*sselam;

for i=1:k
lamdraw(i,1)=gamm_rnd(1,1,.5*vlam1(i,1),.5*sselam(i,1));
lamdraw(i,1)=1/lamdraw(i,1); 
end  

Ht = zeros(k,k+1);
A = sqrt(hdraw)*diag(sqrt(lamdraw));
Ht(:,2:k+1)=A;
Ft=Ht;

    if irep>s0
        %after discarding burnin, store all draws
        h_ = [h_ hdraw];
        lam_ = [lam_ lamdraw]; 
    end
end


alldraws = [h_' lam_'];
%The function momentg is taken from LeSage's toolbox
%it inputs all Gibbs draws and produces posterior
%mean, standard deviation, nse and rne
%it calculates what the book calls S(0) in various ways
%see momentg.m for more details
result = momentg(alldraws);
means=[result.pmean]';
stdevs=[result.pstd]';
nse=[result.nse]';
nse1=[result.nse1]';
nse2=[result.nse2]';
nse3=[result.nse3]';
%calculate Geweke convergence diagnostic based on first .1
%and last .4 of draws
idraw1= round(.1*s1);
result = momentg(alldraws(1:idraw1,:));
meansa=[result.pmean]';
nsea=[result.nse1]';

idraw2= round(.6*s1)+1;
result = momentg(alldraws(idraw2:s1,:));
meansb=[result.pmean]';
nseb=[result.nse1]';

cd = (meansa - meansb)./(nsea+nseb);

'posterior mean, standard deviation and convergence diagnostic, CD'
'h followed by lambda'
[means stdevs cd]

'nse assuming no, .04, .08 and .15 autocovariance estimates'
'h followed by lambda'
[nse nse1 nse2 nse3]

figure(1)
subplot(2,1,1);
hist(lam_(1,:)',20)
title('Posterior Density for \lambda_{0}')
%xlabel('\lambda_{0}')
subplot(2,1,2);
hist(lam_(2,:)',20)
title('Posterior Density for \lambda_{1}')
%xlabel('\lambda_{1}')
%the following material I used for some runs with p=4
%put it back in if you want to experiment with larger values for p
%figure(3)
%hist(lam_(3,:)',20)
%title('Figure 8.2: Posterior Density for \lambda_{2}')
%xlabel('\lambda_{2}')
%figure(4)
%hist(lam_(4,:)',20)
%title('Figure 8.2: Posterior Density for \lambda_{3}')
%xlabel('\lambda_{3}')
%figure(5)
%hist(lam_(5,:)',20)
%title('Figure 8.2: Posterior Density for \lambda_{4}')
%xlabel('\lambda_{4}')