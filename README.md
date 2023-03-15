# ChebfunN
Here is the MATLAB code ChebfunN used in my master thesis, that is a generalization of Chebfun3F for multivariate functions. 
The user can find several function scripts which will all be called in the main chebfunN code. The inputs that must be given to the latter are a
multivariate function f and its number of variables N. For example:

f = @(x,y,z,k) exp(x.*y.*z.*k);

N = 4;

cfN = chebfunN (f,N)

