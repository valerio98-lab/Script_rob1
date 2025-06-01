function S = spline_pchip(t, Q, nGrid)
% Restituisce pos, vel, acc (acc non continua!)
    if nargin < 3, nGrid = 600; end
    [nJ,~] = size(Q);
    tt  = linspace(t(1), t(end), nGrid);

    qq  = zeros(nJ, nGrid);
    vv  = qq;    aa = qq;
    pp  = cell(1,nJ);

    for j = 1:nJ
        pp{j}   = pchip(t, Q(j,:));   % cubic Hermite monotona
        qq(j,:) = ppval(pp{j},            tt);
        vv(j,:) = ppval(fnder(pp{j},1),   tt);   % prima derivata
        aa(j,:) = ppval(fnder(pp{j},2),   tt);   % seconda derivata
    end

    S.pp   = pp;          % cell array dei polinomi
    S.t    = tt;
    S.q    = qq;
    S.v    = vv;
    S.acc  = aa;
end
