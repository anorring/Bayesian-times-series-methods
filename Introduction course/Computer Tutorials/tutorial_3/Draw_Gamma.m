function y = Draw_Gamma(a,b)

% Generate Gamma random draws using the shape-rate parameterization

y = gamrnd(a,1/b);