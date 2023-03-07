function cfN = chebfunN (f, N)
%input: function f and number of variables N
% Initialize

n = linspace(17,17,N);   %coarseResolution
m = n;                   %fineResolution
r = linspace(6,6,N);     %rank
tol = 1e-6;
tol1 = 1e-5;
errmax = Inf;
pref = chebfunpref(); 
cheb = @(i,n) -cos((i-1).*pi/(n-1)); 
reffun = @(n) floor(sqrt(2)^(floor(2*log2(n)) + 1)) + 1;



%% main loop
while errmax > tol1
%% Phase 1
   happyPhase1 = 0;
   while ~happyPhase1 
       
       J = cell(N,1);
       for i=1:N
       J{i} = initializeIndexRandomly(r(i), n(i));
       end
       
       C = cell(N,1);
       for i=1:N
           C{i} = cheb(J{i}, n(i));
       end
       
        
        U = cell(N,1);
        I2 = cell(N,1);
        M = cell(N,1);
        for i=1:N
            Ci = C;
            Ci{i} = -cos(((1:n(i))-1).*pi/(n(i)-1));
            T = constructTensor(f,Ci);
            M{i} = constructMatrix(C,i,n,f);
            [U{i}, ~, ~, I,I2{i}] = ACA(M{i}, tol, n(i)); 
            r(i) = size(I,2);
            J{i} = I;
            C{i} = cheb(J{i}, n(i));
        end
        
            
        refineFlag = 0;
        for i=1:N
            while r(i)*2*sqrt(2) > n(i)
                n(i) = reffun(n(i)); 
                refineFlag = 1;
            end
        end
        if refineFlag == 0
           happyPhase1 = 1;
           break
        end
       
   end

   %% Phase 2

%     for i=1:N
%         if size(J{i},i) == 0 
%                 cfN.U{i} = chebfun(zeros([n(i),1]), pref);
%                 cfN.C = 0;
%         else
        m = n;
        resolved = zeros(N,1); 
        Uf = cell(N,1);
        for l = 1:N
            Uf{l} = U{l};
            ct2 = createCT2(Uf{l});
            resolved(l) = happinessCheck(ct2, [], ct2.coeffs, [], pref);
            if ~resolved(l)
                m(l) = 2*m(l)-1;
            end
        end
            
        % Add function evaluations and check again
        while any(~resolved)
            for p=1:N
                if resolved(p)
                    J{p} = 2*J{p}-1;
                end            
            end
                
            for q=1:N
                if ~resolved(q)
                    M = constructMatrix(C,q,m,f);
                    [~, ~, ~, ~,I2{q}] = ACA(M, tol, m(q));
                    Uf{q} = M(:,I2{q});
                end
                ct2 = createCT2(Uf{q});
                resolved(q) = happinessCheck(ct2, [], ct2.coeffs, [], pref);
                if ~resolved(q)
                    m(q) = 2*m(q)-1;
                end
            end
        end
  
    
            
            
        
    %% Phase 3
        
        % Compute factor matrices
        for k = 1:N
            [Q{k},~] = qr(Uf{k},0);
            K{k} = DEIM(Q{k});    
            U{k} = Q{k}/Q{k}(K{k},:);
        end
        
        % Store factor matrices
         cfN.U = cell(N,1);
         for k = 1:N
             cfN.U{k} = chebfun(U{k}, pref);
         end
                
         % Compute the core
         C = cell(N,1);
         for k=1:N
             C{k} = cheb(K{k}, n(k));       
         end
         cfN.C = constructTensor(f,C);
 
%    end
%   end
%% Calcolo dell'errore di approssimazione
    coeffs = cell(N,1);
    for i= 1:N 
        coeffs{i} = chebcoeffs(full(cfN.U{i}));
    end
    
    A = cfN.C;
    
    for i=1:N
        A = tmprod(A, {coeffs{i}}, i);
    end
    cfN.A = A;
    
    f_approx_vec = zeros(1, length(size(A)));
    fval_vec = zeros(1, length(size(A)));
    err_vec = zeros(1, length(size(A)));
    for i = 1:30
        X = Halton(N,length(size(A)),30);
        f_approx = funapprox(A,X(i,:));
        f_approx_vec(i) = f_approx;
        args = num2cell(X(i,:));
        fval = f(args{:});
        fval_vec(i) = fval;
        err_vec(i) = abs(fval_vec(i) - f_approx_vec(i));
    end
    errmax= max(err_vec);
    
    if any(errmax > tol1)
        for j=1:N
            n(j) = floor(sqrt(2)^(floor(2*log2(n(j))) + 1)) + 1;
            r(j) = 2*r(j);
        end
    elseif any(errmax < tol1)
        break
    end
end
end