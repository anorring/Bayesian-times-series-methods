% Exercise 4
%--------------------------------------------------------------------
% This m-file uses the MATLAB random number generators, and compares
% numerical means and standard deviations to known values. 
%--------------------------------------------------------------------
clear all;
clc;

num_draws = 100; %increase the number of draws to get a better fit to the distribution
%----------------------
%Do results for uniform random number generator
% mean = 0.5, variance = 1/12
%--------------------------
disp('RESULTS - MEAN, STD');
disp('-----------------------------------------------------------------');
disp(' ');
disp('UNIFORM');
tempp = rand(num_draws,1);
unifans = [mean(tempp) std(tempp)];
trueans = [.5 (1/sqrt(12))];
disp([unifans; trueans])
disp(' ');
hist(tempp);
clear tempp;
%----------------------
%Do results for Normal random number generator
%--------------------------
disp('STANDARD NORMAL');
tempp = randn(num_draws,1);  % the drawing command
normans = [mean(tempp) std(tempp)];
trueans = [1e-100 1];
disp([normans; trueans])
disp(' ');
hist(tempp,30);
clear tempp;
%----------------------
%Do results for Student-t random number generator
%--------------------------
disp('STUDENT-T(3)');
tempp = trnd(3,num_draws,1); % the drawing command
tans = [mean(tempp) std(tempp)];
trueans = [1e-100 sqrt(3)];
disp([tans; trueans])
disp(' ');
hist(tempp,30);
clear tempp;
%----------------------
%Do results for Beta random number generator
%--------------------------
disp('BETA(3,2)');
tempp = betarnd(3,2,num_draws,1); % the drawing command
betaans = [mean(tempp) std(tempp)];
trueans = [ (3/(3+2))  ( 1 / 5 )];
disp([betaans; trueans])
disp(' ');
hist(tempp,30);
clear tempp;
%----------------------
%Do results for Exponential random number generator
%--------------------------
disp('Exponential (5)');
tempp = exprnd(5,num_draws,1);
expans = [mean(tempp) std(tempp)];
trueans = [5 5];
disp([expans; trueans])
disp(' ');
hist(tempp);
clear tempp;
%----------------------
%Do results for Chi-Square random number generator
%--------------------------
disp('Chi-Square (3)');
tempp = chi2rnd(3,num_draws,1);
chi2ans = [mean(tempp) std(tempp)];
trueans = [3 sqrt(6)];
disp([chi2ans; trueans])
disp(' ');
hist(tempp,30);
clear tempp;
%----------------------
%Do results for Gamma random number generator
%--------------------------
disp('Gamma(4,2)');
tempp = gamrnd(4,2,num_draws,1);
gamans = [mean(tempp) std(tempp)];
trueans = [8 4];
disp([gamans; trueans])
disp(' ');
hist(tempp,30);
%clear tempp;
