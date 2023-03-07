function ct2 = createCT2(W)
data.hscale = 1;
data.vscale = max(abs(W(:)));
ct2 = chebtech2(W,data);  %Approximate smooth functions on [-1,1] with Chebyshev interpolants
ct2.coeffs = sum(abs(ct2.coeffs), 2);
end