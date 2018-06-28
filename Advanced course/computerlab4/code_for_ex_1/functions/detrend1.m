function [xc,xf] = detrend1(y,q)
% -- Procedure for producing one-sided detrended versions of a 
%    series, using a white-noise + I(2) model.  This produces a 
%    very smooth version of the trend component.  The two-sided 
%    version of the model produces the HP filter.
% 
%    A nice description of what this does is given in 
%    Harvey and Jaeger, JAE, July-Sep. 1993, pp. 231-248
% 
%    Inputs:
%    y = series to be detrended
%    q = relative variance of I(2) component
%        Note: HP quarterly uses q=.000625 (Kydland Prescott)
%              for monthly data a value of q=.00000075
%              matches quarterly gain at 50%, 80% and 90% periods

%  Translated from the Gauss procs of Stock&Watson(2005),'Implications of
%  dynamic factor models for VAR analysis'
%  Dimitris Korobilis, June 2007


vague = 1e+4;  %Vague Prior on I(2) component

%-- Initialize System Matrices -- 
f = zeros(2,2);
f(1,1) = 2;
f(1,2) = -1;
f(2,1) = 1;
x = zeros(2,1);
p = vague*ones(2,2);
p(1,1) = vague+q;
xf = zeros(size(y,1),1);

for i = 1:size(y,1);
    x=f*x;              %x(t/t-1)
    p=f*p*f';
    p(1,1)=p(1,1) + q;  %p(t/t-1)
    if isnan(y(i)) == 0
        h=p(1,1)+1;         %variance of y
        e=y(i)-x(1);        %innovation
        k=p(:,1)/h;         %kalman gain
        x=x+k*e;            %x(t/t)
        p=p-(k*p(1,:));     %p(t/t)
    end
    xf(i)=x(1);
end
xc=y-xf;
