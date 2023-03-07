function indices = DEIM(U)
indices = [];
[~, I] = max(abs(U(:,1)));
indices = [indices,I];
for k = 2:size(U,2)
    c = U(indices,1:(k-1)) \ U(indices,k);
    r = U(:,k) - U(:,1:(k-1))*c;
    [~, I] = max(abs(r));
    indices = [indices,I];
end
end