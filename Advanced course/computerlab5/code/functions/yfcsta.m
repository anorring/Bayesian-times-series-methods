function [yf] = yfcsta(y,tcode,nph)

%  -- Transform series to form series to be forecast 
%      Note input is transformed series instead of raw series
%      This is useful when transformed series contains missing
%      values because outlier ajustments -- 
%      
%      Input:
%      y == transformed series
%      tcode == transformation code
%      nph == forecast horizon
% 
%   -- Tcodes:
%             1 Level
%             2 First Difference
%             3 Second Difference
%             4 Log-Level
%             5 Log-First-Difference
%             6 Log-Second-Difference
% 
%  -- Important Note -- This produces series that is shifted forward nph
%     ahead

n = size(y,1);
yf = zeros(n,1);


% -- Transform and Shift Forward -- 
if tcode==1 || tcode==4   
    yf(1:rows(y)-nph)=y(1+nph:rows(y));
elseif tcode==2 || tcode==5
    for t = 1: n-nph
        yf(t)=sum(y(t+1:t+nph));
    end
elseif tcode==3 || tcode==6
    for t = 1:n-nph
        yf(t)=sum(cumsum(y(t+1:t+nph)));
    end
else
    error('Invalid Transformation Code in yfcst'); 
end

% Divide by nph unless level is being forecast
yf = yf/nph;
if tcode==1 || tcode==4
    yf = yf*nph;
end

yf = yf(1:end-nph,:);