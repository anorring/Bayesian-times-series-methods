function [yf] = yfcst(x,tcode,nph)

% -- Transform series to form series to be forecast -- 
%      
%    Input:
%    x == raw series
%    tcode == transformation code
%    nph == forecast horizon
% 
% -- Tcodes:
%           1 Level
%           2 First Difference
%           3 Second Difference
%           4 Log-Level
%           5 Log-First-Difference
%           6 Log-Second-Difference
%          16 Log(1-L)(1-L^12) forecast deviation
%          17 Log(1-L)(1-L^12) forecast deviation (same as 16)
%            
% -- Important Note -- This produces series that is shifted forward nph ahead

%  Translated from the Gauss procs of Stock&Watson(2005),'Implications of
%  dynamic factor models for VAR analysis'
%  Dimitris Korobilis, June 2007

small=1.0e-06;
n=size(x,1);
yf=zeros(n,1);


% -- Logs or Levels as appropriate -- 
if tcode <= 3
    y=x;
elseif (tcode >= 4) && (tcode <= 6)
    if min(x) < small
        y=NaN;
    end
    y=log(x);
elseif (tcode == 16) || (tcode == 17)
    if min(x) < small
        y=NaN;
    end
    y=log(x);
else
    disp('Invalid Transformation Code in yfcst')
    disp('Tcode = ')
    disp(tcode)
    error('Processing Stops')
end

% -- Transform and Shift Forward -- 
if (tcode == 1) || (tcode == 4)
    yf(1:size(y,1)-nph) = y(1+nph:size(y,1));
elseif (tcode == 2) || (tcode == 5)
    yf(1:size(y,1)-nph) = y(1+nph:size(y,1)) - y(1:size(y,1)-nph);
elseif (tcode == 3) || (tcode == 6)
    yf(2:size(y,1)-nph) = (y(2+nph:size(y,1)) - y(2:size(y,1)-nph)) - nph*(y(2:size(y,1)-nph)-y(1:size(y,1)-nph-1));
elseif (tcode == 16) || (tcode == 17)
    yf(13:size(y,1)-nph) = (y(13+nph:size(y)) - y(13+nph-12:size(y,1)-12)) - (y(13:size(y,1)-nph)-y(13-12:size(y,1)-nph-12)); 
end
