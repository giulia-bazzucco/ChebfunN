function T = constructTensor (f,C)
k = numel(C); 
vec = cell(k,1);
n = cellfun(@numel, C);
for j = 1:k
    linind = 1:n(j);
    vec{j} = C{j}(linind);
end
[vec{:}] = ndgrid(vec{:});
T = f(vec{:});
end