function P = poly434(t, Q, plotFlag)
% poly434  –  traiettoria 4-3-4 (quartica – cubica – quartica)
%            + grafico opzionale integrato
%
% USO -----------------------------------------------------------------
% % dati (2 giunti, 4 nodi)
% t = [0  1.2  2  3];
% Q = [ 0  30  70  40 ;
%       5  20  35  10 ];
% 
% % ► calcola e mostra subito il grafico
% poly434(t, Q);
% 
% % ► solo i coefficienti, nessuna finestra
% P = poly434(t, Q, false);

%   poly434(t, Q)                 % disegna subito la traiettoria
%   P = poly434(t, Q)             % restituisce solo i coefficienti
%   P = poly434(t, Q, false)      % idem, nessun plot
%   poly434(t, Q, true)           % forzatura: calcola + disegna
%
% IN ------------------------------------------------------------------
%   t        1×4     tempi nodi reali [t0 t1 t2 tf]
%   Q        nJ×4    posizioni dei giunti ai nodi
%   plotFlag bool    (default: true se nessun output richiesto)
%
% OUT -----------------------------------------------------------------
%   P (struct) con campi
%       .C1  nJ×5   coeff. quartica:  q_L(τ) = Σ C1(i) τ^(i-1)   τ∈[0,1]
%       .C2  nJ×4   coeff. cubica   :  q_T(τ) = Σ C2(i) τ^(i-1)
%       .C3  nJ×5   coeff. quartica:  q_S(τ) = Σ C3(i) τ^(i-1)

    if nargin < 3, plotFlag = nargout==0; end

    [nJ,N] = size(Q);
    if N~=4, error('Servono esattamente 4 nodi (t,Q).'); end

    t0=t(1); t1=t(2); t2=t(3); tf=t(4);
    h1=t1-t0; h2=t2-t1; h3=tf-t2;

    % ---- velocità costante nel plateau (valida per TUTTI i giunti) ----
    v = (Q(:,3)-Q(:,2))./h2;                           % nJ×1

    % ---- coeff. quartica 1 (lift-off) --------------------------------
    a0 = Q(:,1);                           a1 = 0;
    a2 = 0;
    a3 =  4*(Q(:,2)-Q(:,1)) - v*h1;
    a4 = -3*(Q(:,2)-Q(:,1)) + v*h1;
    C1 = [a0 a1*ones(nJ,1) a2*ones(nJ,1) a3 a4];

    % ---- coeff. cubica (travel) --------------------------------------
    b0 = Q(:,2);
    b1 = v*h2;
    S1 = Q(:,3)-Q(:,2) - b1;
    S2 = v*h2 - b1;
    b3 = S2 - 2*S1;
    b2 = S1 - b3;
    C2 = [b0 b1 b2 b3];

    % ---- coeff. quartica 2 (set-down) --------------------------------
    c0 = Q(:,3);
    c1 = v*h3;
    D  = Q(:,4)-Q(:,3) - c1;
    c4 =  3*D + 2*c1;
    c3 = (c1 - 8*c4)/3;
    c2 = -3*c3 - 6*c4;
    C3 = [c0 c1 c2 c3 c4];

    % ---- pacchetto di uscita ----------------------------------------
    P.C1 = C1;  P.C2 = C2;  P.C3 = C3;

    % ---- grafico opzionale ------------------------------------------
    if plotFlag
        nGrid = 600;
        tt = linspace(t0, tf, nGrid);
        qq = zeros(nJ, numel(tt));
        for k = 1:numel(tt)
            ti = tt(k);
            if ti <= t1             % tratto lift-off
                tau = (ti-t0)/h1;
                for j=1:nJ
                    qq(j,k)=polyval(flip(C1(j,:)), tau);
                end
            elseif ti <= t2         % tratto travel
                tau = (ti-t1)/h2;
                for j=1:nJ
                    qq(j,k)=polyval(flip(C2(j,:)), tau);
                end
            else                    % tratto set-down
                tau = (ti-t2)/h3;
                for j=1:nJ
                    qq(j,k)=polyval(flip(C3(j,:)), tau);
                end
            end
        end

        col = lines(nJ);
        figure('Name','4-3-4 trajectory'), clf, hold on, grid on
        for j = 1:nJ
            plot(tt, qq(j,:), 'Color', col(j,:), 'LineWidth',1.6, ...
                 'DisplayName', sprintf('joint %d',j));
            plot(t, Q(j,:), 'o', 'Color', col(j,:), ...
                 'MarkerFaceColor','w', 'HandleVisibility','off');
        end
        xlabel('t'), ylabel('q_j'),
        title('Traiettoria 4-3-4   (quartica – cubica – quartica)')
        legend('Location','best')
    end
end
