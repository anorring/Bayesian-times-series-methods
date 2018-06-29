
% This procedure follows De Jong and Shephard (1995), "The Simulation
% Smoother for Time Series Models", Biometrika 82, 339-350. Ft is assumed
% to matrix of full row rank. We use equation (3). Output of the procedure
% is a draw eta(t)
%Note that I am using notation from the DeJong and Shephard paper which differs a bit
%from that used in the book
%This code does not allow for regression effects (i.e. no X)
%There are also some things set up for univariate case, alterations required for multivariate
%
%applies the Kalman filter to all t

function nt = djs(yt,Zt,Gtt,Tt,Ht,Ft,sig);

T  = size(yt,1);
N  = size(Tt,1);


% First, run Kalman Filter 

at = zeros(N,1);
Pt=Ht*Ht';
et = zeros(T,1);
DIt= zeros(T,1);
Ltt= zeros(T,N*N);
rowh=size(Ht,1);
colh=size(Ht,2);
Jtt= zeros(T,N*colh);

for i=1:T
Ztt=Zt(i,:);
    et(i,1) = yt(i,1)-Ztt*at;
    D       = Ztt*Pt*Ztt'+Gtt*Gtt';
    DIt(i,1)  = 1/D;
    Kt      = (Tt*Pt*Ztt'+Ht*Gtt')*DIt(i,1);
    L       = Tt-Kt*Ztt;
    Jt      = Ht-Kt*Gtt;
    at      = Tt*at+Kt*et(i,1);
    Pt      = Tt*Pt*L'+Ht*Jt';
    Ltt(i,:)= reshape(L,N^2,1)';
    Jtt(i,:)= reshape(Jt,rowh*colh,1)';
   
end

% Next, use (4) to sample F(t)u(t) 

rt = zeros(N,1);
Ut = zeros(N,N);

rowf=size(Ft,1);
nt = zeros(T,rowf);

i = T;

while i>0
Ztt=Zt(i,:);
    L       = reshape(Ltt(i,:),N,N);
    Jt      = reshape(Jtt(i,:),N,colh);
    Ct      = Ft*Ft'-Ft*Gtt'*DIt(i,1)*Gtt*Ft'-Ft*Jt'*Ut*Jt*Ft';
    CIt     = inv(Ct);
    eps     = sig*chol(Ct)'*randn(rowf,1);
    nt(i,:) = (Ft*Gtt'*DIt(i,1)*(et(i,:)')+Ft*Jt'*rt+eps)';
    Vt      = Ft*Gtt'*DIt(i,1)*Ztt+Ft*Jt'*Ut*L;
    rt      = Ztt'*DIt(i,1)*et(i,1)+L'*rt-Vt'*CIt*eps;
    Ut      = Ztt'*DIt(i,1)*Ztt+L'*Ut*L+Vt'*CIt*Vt;
    i=i-1;
end

% Finally, in notation of DJS, draw nt(0) 

Ct  = Ft*Ft'-Ft*Ht'*Ut*Ht*Ft';
nt0 = (Ft*Ht'*rt+sig*chol(Ct)'*randn(rowf,1))';


% Note that the algorithm of De Jong/Shephard draws nt(0),...,nt(T).
% To compute the relevant states, only nt(0),...,nt(T-1) are needed.

nt  = [nt0; nt];
