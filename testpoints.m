function err = testpoints(N,A,f)
%input: number of variables N, tensor A and function f
f_approx_vec = zeros(1, length(size(A)));
fval_vec = zeros(1, length(size(A)));
err_vec = zeros(1, length(size(A)));

for i = 1:100
    X = Halton(N,length(size(A)),100);
    f_approx = funapprox(A,X(i,:));
    f_approx_vec(i) = f_approx;
    args = num2cell(X(i,:));
    fval = f(args{:}); 
    fval_vec(i) = fval;
    err_vec(i) = abs(fval_vec(i) - f_approx_vec(i));
end
err = norm(err_vec, "fro"); 
end