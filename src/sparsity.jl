function  expandLambda( lambda, nf )

    # expandLambda  expand a parent array
    # expandLambda(lambda,nf) calculates the expanded parent array, for use in
    # sparse factorization algorithms, from a given parent array and an array
    # of joint motion freedoms.  nf(i) is the degree of motion freedom allowed
    # by joint i.
    
    N = length(lambda)
    n = sum(nf)
    
    newLambda = collect(0:(n-1))
    
    map(1) = 0;				# (matlab won't let us use map(0))
    for i = 1:N
      map(i+1) = map(i) + nf(i);
    end
    
    for i = 1:N
      newLambda[map[i]+1] = map[lambda[i]+1]
    end
    return newLambda
end

    function  [L,D] = LTDL( H, lambda )
    
    # LTDL  factorize H -> L'*D*L exploiting branch-induced sparsity
    # [L,D]=LTDL(H,lambda) returns a unit-lower-triangular matrix L and a
    # diagonal matrix D that satisfy L'*D*L = H, where H is a symmetric,
    # positive-definite matrix having the property that the nonzero elements on
    # row i below the main diagonal appear only in columns lambda(i),
    # lambda(lambda(i)), and so on.  This is the pattern of branch-induced
    # sparsity; and H and lambda can be regarded as the joint-space inertia
    # matrix and parent array of a kinematic tree.  lambda must satisfy
    # 0<=lambda(i)<i for all i.
    
    n = size(H,1);
    
    for k = n : -1 : 1
      i = lambda(k);
      while i ≠ 0
        a = H(k,i) / H(k,k);
        j = i;
        while j ≠ 0
          H(i,j) = H(i,j) - a * H(k,j);
          j = lambda(j);
        end
        H(k,i) = a;
        i = lambda(i);
      end
    end
    
    D = diag(diag(H));
    L = eye(n);
    for i = 2:n
      L[i,1:(i-1)] = H[i,1:(i-1)]
    end
end

    function  L = LTL( H, lambda )
    
    # LTL  factorize H -> L'*L exploiting branch-induced sparsity
    # LTL(H,lambda) returns a lower-triangular matrix L satisfying L'*L = H,
    # where H is a symmetric, positive-definite matrix having the property that
    # the nonzero elements on row i below the main diagonal appear only in
    # columns lambda(i), lambda(lambda(i)), and so on.  This is the pattern of
    # branch-induced sparsity; and H and lambda can be regarded as the
    # joint-space inertia matrix and parent array of a kinematic tree.  lambda
    # must satisfy 0<=lambda(i)<i for all i.
    
    n = size(H,1);
    
    for k = n : -1 : 1
      H(k,k) = sqrt(H(k,k));
      i = lambda(k);
      while i ≠ 0
        H(k,i) = H(k,i) / H(k,k);
        i = lambda(i);
      end
      i = lambda(k);
      while i ≠ 0
        j = i;
        while j ≠ 0
          H(i,j) = H(i,j) - H(k,i) * H(k,j);
          j = lambda(j);
        end
        i = lambda(i);
      end
    end
    
    L = H;
    for i = 1:n-1  				# zero the upper triangle
      L(i,i+1:n) = zeros(1,n-i);
    end
end

    function  y = mpyH( H, lambda, x )
    
    # mpyH  calculate H*x exploiting branch-induced sparsity in H
    # mpyH(H,lambda,x) computes H*x where x is a vector and H is a symmetric,
    # positive-definite matrix having the property that the nonzero elements on
    # row i below the main diagonal appear only in columns lambda(i),
    # lambda(lambda(i)), and so on.  This is the pattern of branch-induced
    # sparsity; and H and lambda can be regarded as the joint-space inertia
    # matrix and parent array of a kinematic tree.  lambda must satisfy
    # 0<=lambda(i)<i for all i.
    
    n = size(H,1);
    y = x;           # to give y same dimensions as x
    
    for i = 1:n
      y(i) = H(i,i) * x(i);
    end
    for i = n:-1:1
      j = lambda(i);
      while j ≠ 0
        y(i) = y(i) + H(i,j) * x(j);
        y(j) = y(j) + H(i,j) * x(i);
        j = lambda(j);
      end
    end
end
    function  y = mpyL( L, lambda, x )
    
    # mpyL  multiply vector by L factor from LTL or LTDL
    # mpyL(L,lambda,x) computes L*x where L is the lower-triangular matrix from
    # either LTL or LTDL and lambda is the parent array describing the sparsity
    # pattern in L.
    
    n = size(L,1);
    y = x;           # to give y same dimensions as x
    
    for i = 1:n
      y(i) = L(i,i) * x(i);
      j = lambda(i);
      while j ≠ 0
        y(i) = y(i) + L(i,j) * x(j);
        j = lambda(j);
      end
    end
end
    function  y = mpyLi( L, lambda, x )
    
    # mpyLi  multiply vector by inverse of L factor from LTL or LTDL
    # mpyLi(L,lambda,x) computes inv(L)*x where L is the lower-triangular
    # matrix from either LTL or LTDL and lambda is the parent array describing
    # the sparsity pattern in L.
    
    n = size(L,1);
    
    for i = 1:n
      j = lambda(i);
      while j ≠ 0
        x(i) = x(i) - L(i,j) * x(j);
        j = lambda(j);
      end
      x(i) = x(i) / L(i,i);
    end
    
    y = x;
end
    function  y = mpyLit( L, lambda, x )
    
    # mpyLit  multiply vector by inverse transpose of L factor from LTL or LTDL
    # mpyLit(L,lambda,x) computes inv(L')*x where L is the lower-triangular
    # matrix from either LTL or LTDL and lambda is the parent array describing
    # the sparsity pattern in L.
    
    n = size(L,1);
    
    for i = n:-1:1
      x(i) = x(i) / L(i,i);
      j = lambda(i);
      while j ≠ 0
        x(j) = x(j) - L(i,j) * x(i);
        j = lambda(j);
      end
    end
    
    y = x;
end
    function  y = mpyLt( L, lambda, x )
    
    # mpyLt  multiply vector by transpose of L factor from LTL or LTDL
    # mpyLt(L,lambda,x) computes L'*x where L is the lower-triangular matrix
    # from either LTL or LTDL and lambda is the parent array describing the
    # sparsity pattern in L.
    
    n = size(L,1);
    y = x;           # to give y same dimensions as x
    
    for i = 1:n
      y(i) = L(i,i) * x(i);
      j = lambda(i);
      while j ≠ 0
        y(j) = y(j) + L(i,j) * x(i);
        j = lambda(j);
      end
    end
end