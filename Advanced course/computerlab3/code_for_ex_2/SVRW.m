function h = SVRW(ystar,h,h0,sigh2)
T = length(h);
    % normal mixture
pj = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mj = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819]...
    - 1.2704;  % caution: means are adjusted!
sigj2 = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];
sigj = sqrt(sigj2);

    % sample S from a 7-point distrete distribution
temprand = rand(T,1);
q = repmat(pj,T,1).*normpdf(repmat(ystar,1,7),...
    repmat(h,1,7)+repmat(mj,T,1),repmat(sigj,T,1));
q = q./repmat(sum(q,2),1,7);
S = 7 - sum(repmat(temprand,1,7)<cumsum(q,2),2) + 1;

    % sample h
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
HH = H'*H;
d_s = mj(S)'; 
iSig_s = sparse(1:T,1:T,1./sigj2(S));
Kh = HH/sigh2 + iSig_s;
h_hat = Kh\(h0/sigh2*HH*ones(T,1) + iSig_s*(ystar-d_s));
h = h_hat + chol(Kh,'lower')'\randn(T,1);
end