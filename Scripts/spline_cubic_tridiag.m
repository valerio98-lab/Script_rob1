function S = spline_cubic_tridiag(t, Q, v1, vN)
% spline_cubic_tridiag
%  Cubic C² spline per un robot a nJ giunti e N nodi.
%  - risolve il sistema tridiagonale in blocco per tutti gli assi
%  - restituisce velocità, accelerazioni e coefficienti per tratto
%
% USO --------------------------------------------------------------
%   S = spline_cubic_tridiag(t, Q)            % v1 = vN = 0   (default)
%   S = spline_cubic_tridiag(t, Q, v1, vN)    % v1,vN vettori nJ×1

% tempi dei 4 nodi
% t = [1 2 2.5 4];
% 
% % tre giunti, quattro nodi (3 × 4)
% Q = [45  90  -45   45  ;    % joint 1
%       0  20   10   10  ;    % joint 2
%       5   5    8   12 ];    % joint 3
% 
% S = spline_cubic_tridiag(t, Q);    % v1 = vN = 0 per tutti
% 
% % velocità interne stampate
% disp(S.v)
% 
% % valutazione del giunto 2 in t = 2.3 s
% k = find(t<=2.3,1,'last');               % tratto corrente
% tau = 2.3 - t(k);
% coeff = squeeze(S.C(2,k,:));             % coeffic. joint #2 tratto k
% q23 = polyval(flipud(coeff), tau);
% fprintf('q2(2.3 s) = %.3f\n', q23);

%
% IN ---------------------------------------------------------------
%   t : 1×N   tempi (strettamente crescenti)
%   Q : nJ×N  ogni riga = valori di un giunto ai N nodi
%   v1,vN (fac.) : velocità iniziale / finale (nJ×1)
%
% OUT  S  (struct)
%   .v    nJ×N           velocità ai nodi
%   .acc  nJ×N           accelerazioni ai nodi
%   .C    nJ×(N-1)×4     coefficienti cubici [a0 a1 a2 a3]

    if nargin < 3, v1 = []; end
    if nargin < 4, vN = []; end

    [nJ, N] = size(Q);
    t = t(:).';                           % riga
    if numel(t) ~= N                     % check rapido
        error('Dimensioni non coerenti fra t e Q');
    end

    h = diff(t);                         % 1×(N-1)

    % ---------- matrice tridiagonale A (N-2 × N-2) ----------
    main = 2*(h(1:end-1)+h(2:end));
    A = diag(main) + diag(h(2:end-1),1) + diag(h(2:end-1),-1);

    % ---------- termine noto RHS (N-2 × nJ) -----------------
    dq = diff(Q,1,2);                    % nJ×(N-1)
    R1 = h(1:end-1) .* (dq(:,2:end) ./ h(2:end));
    R2 = h(2:end)   .* (dq(:,1:end-1) ./ h(1:end-1));
    RHS = 3 * (R1 + R2).';               % trasposta → (N-2)×nJ

    % correzione se v1/vN fissate (clamped)
    if ~isempty(v1)
        RHS(1,:)  = RHS(1,:)  - h(1) * v1(:).';
    end
    if ~isempty(vN)
        RHS(end,:) = RHS(end,:) - h(end) * vN(:).';
    end

    % ---------- soluzione tridiagonale simultanea -----------
    Vint = A \ RHS;                      % (N-2)×nJ

    % velocità complete
    v = [zeros(1,nJ); Vint; zeros(1,nJ)];   % default v1=vN=0
    if ~isempty(v1), v(1,:)   = v1(:).'; end
    if ~isempty(vN), v(end,:) = vN(:).'; end
    v = v.';                               % nJ×N

    % ---------- coefficienti cubici per ogni giunto ----------
    C = zeros(nJ, N-1, 4);                 % pre-alloc
    for k = 1:N-1
        hk = h(k);
        dk  = (Q(:,k+1)-Q(:,k)) / hk;      % nJ×1
        a0 = Q(:,k);                       % nJ×1
        a1 = v(:,k);
        a2 = (3*dk - 2*v(:,k) - v(:,k+1)) / hk;
        a3 = (2*(Q(:,k) - Q(:,k+1)) + hk*(v(:,k)+v(:,k+1))) / hk^2;
        C(:,k,1) = a0;  C(:,k,2) = a1;
        C(:,k,3) = a2;  C(:,k,4) = a3;
    end

    % ---------- accelerazione nei nodi -----------------------
    acc = zeros(nJ, N);
    acc(:,1) = 2*C(:,1,3);
    for k = 2:N-1
        acc(:,k) = 2*C(:,k,3);
    end
    acc(:,N) = 2*C(:,N-1,3) + 6*C(:,N-1,4).*h(end);

    % ---------- struct di uscita -----------------------------
    S.v   = v;
    S.acc = acc;
    S.C   = C;
end
