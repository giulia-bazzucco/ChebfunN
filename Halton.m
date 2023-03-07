function X = Halton (N,d,numpnt)
p = haltonset(d,'Skip',1); 
H = p(1:numpnt,:);
X = zeros(size(H));
for i=1:N
    X(:,i) = H(:,i);
    X(:,i) = 2*X(:,i) -1;
end
end
