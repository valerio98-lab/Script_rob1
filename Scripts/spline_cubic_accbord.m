function S = spline_cubic_accbord(t, Q, acc0, accN)
% spline_cubic_accbord
%   Cubic C² con accelerazione imposta in t(1) e/o t(end)
%   via ghost-knots automatici.

%% UTILIZZO
% t  = [1 2 2.5 4];                 % nodi reali (4)
% Q  = [45 90 -45 45 ;              % giunto 1
%        0 20  10 10];              % giunto 2   -> 2×4
% 
% acc0 = [0 ; 0];                   % voglio acc = 0 al primo nodo
% accN = [0 ; 0];                   % e acc = 0 all’ultimo nodo
% 
% S = spline_cubic_accbord(t, Q, acc0, accN);
% 
% % velocità calcolate
% disp('v ai nodi reali:')
% disp(S.v)
% 
% % accelerazione ai nodi (vedrai 0 a inizio/fine)
% disp('acc ai nodi reali:')
% disp(S.acc)
%%
%
% IN -----------------------------------------------------------------
%   t     1×N     tempi nodi reali (strettamente crescenti)
%   Q     nJ×N    posizioni giunto ai nodi reali
%   acc0  nJ×1    accelerazione desiderata in t(1)    (facoltativo)
%   accN  nJ×1    accelerazione desiderata in t(end)  (facoltativo)
%
% OUT ----------------------------------------------------------------
%   S  = struct identico a spline_cubic_tridiag
%        (S.v, S.acc, S.C)   ma riferito **solo** all’intervallo reale
%
% Dipende da spline_cubic_tridiag  (stesso file o path MATLAB)

    if nargin < 3, acc0 = []; end
    if nargin < 4, accN = []; end

    [nJ, N] = size(Q);
    t = t(:).';                              % riga

    % ---- costruiamo vettori estesi ---------------------------------
    tExt = t;
    QExt = Q;

    % --- nodo fantasma PRIMA di t1 ----------------------------------
    if ~isempty(acc0)
        h0 = t(2) - t(1);                    % stessa ampiezza del 1° intervallo
        q0 = Q(:,1) + h0^2 * acc0 + 2*Q(:,1) - Q(:,2);  % q0 = h0^2*acc + 2q1 - q2
        tExt = [t(1)-h0 , tExt];
        QExt = [q0     ,  QExt];
    end

    % --- nodo fantasma DOPO tN --------------------------------------
    if ~isempty(accN)
        hN = t(end) - t(end-1);
        qG = Q(:,end) + hN^2*accN + 2*Q(:,end) - Q(:,end-1);
        tExt = [tExt , t(end)+hN];
        QExt = [QExt , qG];
    end

    % ---- spline cubica standard sui nodi estesi --------------------
    Sfull = spline_cubic_tridiag(tExt, QExt);   % v1 = vN = 0 di default

    % ---- ritagliamo le parti che interessano l’intervallo reale ----
    if ~isempty(acc0)
        idx0 = 2;          % si parte dal 2° nodo dell'estesa
    else
        idx0 = 1;
    end
    if ~isempty(accN)
        idxN = size(tExt,2)-1;  % penultimo
    else
        idxN = size(tExt,2);
    end

    S.v   = Sfull.v(:, idx0:idxN);
    S.acc = Sfull.acc(:, idx0:idxN);
    S.C   = Sfull.C(:, idx0-1:idxN-1, :);
end
