function X_out = SURform(X)
[r,c] = size( X );
idi = kron((1:r)',ones(c,1));
idj = (1:r*c)';
X_out = sparse(idi,idj,reshape(X',r*c,1));
end