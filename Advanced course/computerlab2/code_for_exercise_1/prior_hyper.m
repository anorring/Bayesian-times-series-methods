% Define hyperparameters for SSVS prior


    n = K*M; % Total number of parameters (size of vector alpha)
    % mean of alpha
    a_prior = zeros(n,1);
    
    % SSVS variances for alpha
    tau_0=0.01*ones(n,1);
    tau_1=100*ones(n,1);
    
    % Priors on SIGMA
  
        % SSVS variances for non-diagonal elements of SIGMA     
        kappa_0 = 0.01; % Set kappa_[0ij], kappa_[1ij]
        kappa_1 = 100;
   
        % Hyperparameters for diagonal elements of SIGMA (Gamma density)
        a_i = .01;
        b_i = .01;
    
    
    % Hyperparameters for Gamma ~ BERNOULLI(m,p_i), see eq. (14)
    p_i = .5;
    
    % Hyperparameters for Omega_[j] ~ BERNOULLI(j,q_ij), see eq. (17)
    q_ij = .5;
    
    % Initialize Gamma and Omega vectors
    gammas = ones(n,1);       % vector of Gamma
    omega=cell(1,M-1);
    for kk_1 = 1:(M-1)
        omega{kk_1} = ones(kk_1,1);	% Omega_j
    end
    
    % Set space in memory for some vectors that we are using in SSVS
    gamma_draws = zeros(nsave,n); % vector of gamma draws
    omega_draws = zeros(nsave,.5*M*(M-1)); %vector of omega draws    
