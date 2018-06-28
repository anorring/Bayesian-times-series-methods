%--------------------------------------------------------------------------
% Bayesian VAR estimation using the Minnesota Prior
%The VAR model is:
%
%     Y(t) = A0 + Y(t-1) x A1 + ... + Y(t-p) x Ap + e(t)
%
% The Model
%
%                   Y(t) = X(t) x A + e(t)
%
% where e(t) ~ N(0,SIGMA), and A summarizes all parameters. Note that we
% also use the vector a which is defined as a=vec(A). 
%
% Code written by Dimitris Korobilis and modified by Gary Koop
%--------------------------------------------------------------------------

clear all;

%------------------------------LOAD DATA-----------------------------------
% Load Quarterly US data on inflation, unemployment and interest rate, 

[Yraw,Dates] = xlsread('US_macrodata.csv');


% If you want to use this code with your data set, name the data you load 'Yraw'
%Note that 'Yraw' is a matrix with T rows by M columns,
% where T is the number of time series observations and M is the number of VAR dependent variables.
%----------------------------PRELIMINARIES---------------------------------
% Define specification of the VAR model
constant = 1;        % 1: if you desire intercepts, 0: otherwise 
p = 1;               % Number of lags on dependent variables


%--------------------------DATA HANDLING-----------------------------------
% Get initial dimensions of dependent variables
[Traw M] = size(Yraw);

   Y1 = Yraw;

        
% Generate lagged Y matrix. This will be part of the X matrix
Ylag = mlag2(Y1,p); % Y is [T x M]. ylag is [T x (Mp)]

% Now define matrix X which has all the R.H.S. variables (constant, lags of
% the dependent variable)

if constant
    X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:)];
else
    X1 = Ylag(p+1:Traw,:);  %
end

% Get size of final matrix X
[Traw3 K] = size(X1);


% Form Y matrix accordingly
% Delete first "LAGS" rows to match the dimensions of X matrix
Y1 = Y1(p+1:Traw,:); % This is the final Y matrix used for the VAR

% Traw was the dimension of the initial data. T is the number of actual 
% time series observations of Y and X
T = Traw - p;


    Y = Y1;
    X = X1;


%--------------------------------PRIORS------------------------------------
% First get Ordinary Least Squares (OLS) estimators 
% These are to be used in posterior formulae
A_OLS = inv(X'*X)*(X'*Y); % This is the matrix of regression coefficients
a_OLS = A_OLS(:);         % This is the vector of coefficients, i.e. it holds vec(A_OLS)
SSE = (Y - X*A_OLS)'*(Y - X*A_OLS);
SIGMA_OLS = SSE./(T-K);

%-----------------Prior hyperparameters for bvar model
% Prior mean
    A_prior = 0*ones(K,M);   
    a_prior = A_prior(:);
    
    % Hyperparameters on the Minnesota variance of alpha
    a_bar_1 = 0.5;
    a_bar_2 = 0.25;
    a_bar_3 = 10^2;
    
    % Now get residual variances of univariate p-lag autoregressions.
    %These are used in the Minnesota prior covariance matrix
    
    sigma_sq = zeros(M,1); % vector to store residual variances
    for i = 1:M
        % Create lags of dependent variable in i-th equation
        Ylag_i = mlag2(Yraw(:,i),p);
        Ylag_i = [ones(T,1), Ylag_i(p+1:Traw,:)];
        % Dependent variable in i-th equation
        Y_i = Yraw(p+1:Traw,i);
        % OLS estimates of i-th equation
        alpha_i = inv(Ylag_i'*Ylag_i)*(Ylag_i'*Y_i);
        sigma_sq(i,1) = (1./(T-p+1))*(Y_i - Ylag_i*alpha_i)'*(Y_i - Ylag_i*alpha_i);
    end
    
   
    % Create an array of dimensions K x M, which will contain the K diagonal
    % elements of the prior covariance matrix, in each of the M equations.
    V_i = zeros(K,M);
    % index in each equation which are the own lags
    ind = zeros(M,p);
    for i=1:M
        ind(i,:) = constant+i:M:K;
    end
    for i = 1:M  % for each i-th equation
        for j = 1:K   % for each j-th RHS variable
            if constant==1
                if j==1 % if there is constant, use this code
                    V_i(j,i) = a_bar_3*sigma_sq(i,1); % variance on constant                
                elseif find(j==ind(i,:))>0
                    V_i(j,i) = a_bar_1./(ceil((j-1)/M)^2); % variance on own lags           
                    % Note: the "ceil((j-1)/M)" command finds the associated lag 
                    % number for each parameter
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;                   
                        end
                    end                 % variance on other lags  
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((ceil((j-1)/M)^2)*sigma_sq(ll,1));           
                end
            else   % if no constant is defined, then use this code
                if find(j==ind(i,:))>0
                    V_i(j,i) = a_bar_1./(ceil(j/M)^2); % variance on own lags
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;
                        end                        
                    end                 % variance on other lags  
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((ceil(j/M)^2)*sigma_sq(ll,1));            
                end
            end
        end
    end
 
    
    % Now V is a diagonal matrix with diagonal elements the V_i
    V_prior = diag(V_i(:));  % this is the prior variance of the vector a  
    
   

    
%============================ POSTERIORS ==================================
%==========================================================================
    

%--------- Posterior hyperparameters of ALPHA and SIGMA with Minnesota Prior
    % ******Get all the required quantities for the posteriors       
    V_post = inv( inv(V_prior) + kron(inv(SIGMA_OLS),X'*X) );
    a_post = V_post * ( inv(V_prior)*a_prior + kron(inv(SIGMA_OLS),X'*X)*a_OLS );
    A_post = reshape(a_post,K,M);
     
    

% Print some results

disp('The posterior mean of alpha')
A_post

disp('The posterior st dev of alpha')
reshape((diag(sqrt(V_post))),K,M)

