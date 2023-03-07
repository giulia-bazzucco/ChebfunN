function[Mc, Mr, Mt, rInd, cInd] = ACA (M, tol, maxIter)
Mc = [];
Mr = [];
Mt = [];
rInd = [];
cInd = [];
Moriginal = M;

for iter=1:maxIter
    
    [err,I2] = max(abs(M(:)));
    if isempty(err) || err < tol
        Mc = Moriginal(:,cInd);
        Mr = Moriginal(rInd,:)';
        Mt = Moriginal(rInd,cInd);
        return
    end
    
    [I,J] = ind2sub(size(M), I2); 
    rInd = [rInd, I];
    cInd = [cInd, J];
    

    M = M-M(:,J)*M(I,:)./M(I,J);
end
Mc = Moriginal(:,cInd);
Mr = Moriginal(rInd,:)';
Mt = Moriginal(rInd,cInd);
end