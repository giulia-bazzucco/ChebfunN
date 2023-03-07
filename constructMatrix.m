function M = constructMatrix (C,k,n,f)

Cp = C;
cheb = @(i,n) -cos((i-1).*pi/(n-1));
Cp{k} = cheb(1:n(k), n(k));
N = length(C);

T = constructTensor (f,Cp);
T = permute(T , [k, 1:k-1, k+1:N]); 
M = reshape(T, n(k), []); 

end
