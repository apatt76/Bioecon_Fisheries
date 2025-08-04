function F = computeF(X, p)
    if size(X,2) > 1
        X = X(:,1);  % only use first column if multiple passed
    end

    spe = length(p.nicheweb);
    B = X(1:spe);
    B(B < p.Bext) = 0;

    basal = (sum(p.nicheweb, 2) == 0);
    nonbasal = ~basal;

    w = zeros(spe,1);
    w(nonbasal) = 1 ./ sum(p.nicheweb(nonbasal, B~=0), 2);
    w(w==Inf) = 0;
    w = w * ones(1, spe);
    w = w .* p.nicheweb;


    wBh = w .* (ones(spe,1)*B').^p.h;
    Bsh = (ones(spe,1) * p.Bs').^p.h;
    cBBsh = (p.c .* B * ones(1, spe)) .* Bsh;
    sumwBh = sum(wBh, 2) * ones(1, spe);
    F = wBh ./ (Bsh + cBBsh + sumwBh);


    
end
