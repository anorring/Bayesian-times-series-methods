% function [By,Bx,u,omega,xx]=estvar(ydata,lags,xdata)
% estimates a VAR, parts of this code use Sims's rfvar3.m
% By is equ*var*lags, Bx is equ*nx, omega is ny*ny, xx is (ny*lags+nx)*(ny*lags+nx)
function [By,Bx,u,omega,xx]=estvar(ydata,lags,xdata)

[t,ny]=size(ydata);
nox=isempty(xdata);
if ~nox
   [t1,nx]=size(xdata);
else
   t1=t; nx=0;
end
if t1~=t, disp('Mismatch of x and y data lengths'), end

X=zeros(t,ny,lags);
for i=1:lags
   X(lags+1:t,:,i)=ydata(lags+1-i:t-i,:);
end
X=[X(:,:) xdata];
X=X(lags+1:t,:);
y=ydata(lags+1:t,:);

k=ny*lags+nx;	% # of coefficients
[vl,d,vr]=svd(X,0);	% singular value decomposition, X is t*k
d=1./diag(d);
B=vl'*y;
B=(vr.*repmat(d',k,1))*B;	% this is OLS equ by equ, B is k*ny
u=y-X*B;
omega=u'*u/(t-lags);	% covarinace matrix
xx=vr.*repmat(d',k,1);
xx=xx*xx';
By=B(1:ny*lags,:);
By=reshape(By,ny,lags,ny); % variables, lags, equations
By=permute(By,[3,1,2]); % equations, variables, lags

if nox
   Bx=[];
else
   Bx=B(ny*lags+(1:nx),:)'; % equations*nx
end
