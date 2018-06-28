function yIR = construct_IR(beta,Sig,n_hz,shock)
n = size(Sig,1);
p = (size(beta,1)/n-1)/n;
CSig = chol(Sig,'lower');
tmpZ1 = zeros(p,n); tmpZ = zeros(p,n);
Yt1 = CSig*shock; Yt = zeros(n,1);
yIR = zeros(n_hz,n); yIR(1,:) = Yt1';
for t = 2:n_hz    
        % update the regressors    
    tmpZ = [Yt'; tmpZ(1:end-1,:)];
    tmpZ1 = [Yt1'; tmpZ1(1:end-1,:)];    
        % evolution of variables if a shock hits
    e = CSig*randn(n,1);
    Z1 = reshape(tmpZ1',1,n*p);
    Xt1 = kron(speye(n), [1 Z1]); 
    Yt1 = Xt1*beta + e;    
        % evolution of variables if no shocks hit
    Z = reshape(tmpZ',1,n*p);
    Xt = kron(speye(n), [1 Z]);
    Yt = Xt*beta + e;    
        % the IR is the difference of the two scenarios
    yIR(t,:) = (Yt1-Yt)';     
end    
end