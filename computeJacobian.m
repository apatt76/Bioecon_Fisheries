function J = computeJacobian(fun, X, t, params)

    if size(X,2) > 1
        X = X(:,1);  % only use first column if multiple passed
    end
    
    n = length(X);
    fx = fun(t, X, params);
    delta = 1e-6;
    J = zeros(n);
    for i = 1:n
        X_eps = X;
        X_eps(i) = X_eps(i) + delta;
        fx_eps = fun(t, X_eps, params);
        J(:,i) = (fx_eps - fx) / delta;
    end
end