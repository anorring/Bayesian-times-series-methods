function [Yraw] = transform_hstep(ydata,tcode,nfore)
%TRANSFORM Transform large dataset to stationarity
%This code corrects the number of observations lost from transformations
Yraw = 0*ydata;
for i=1:size(ydata,2)
   Yraw(:,i) = yfcst(ydata(:,i),tcode(i),nfore); 
end
Yraw = Yraw(1:end-nfore,:);
end

