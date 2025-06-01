function S = spline_plot_and_compute(t, Q, method, varargin)
% % La funzione calcola direttamente la spline ma in più agli altri script
% la plotta anche. Quando avrò il tempo li renderò una cosa unica

% method: 'cubicC2' (default) | 'pchip'

%% UTILIZZO 

% % dati di esempio (3 giunti, 4 nodi)
% t = [1  2  2.5  4];
% Q = [ 40  70 -50  40 ;
%        0  10  20  10 ;
%        5   5   8  12 ];
% 
% % ---- Caso 1: v1 = vN = 0 (cubica "standard")
% S = spline_plot_and_compute(t, Q);
% 
% % ---- Caso 2: impongo accelerazione iniziale e finale = 0
% acc0 = zeros(size(Q,1),1);     % [0;0;0]
% accN = zeros(size(Q,1),1);     % idem
% S = spline_plot_and_compute(t, Q, acc0, accN);

%%
% quick_spline_demo(t,Q,...) –
%    Se passi acc0/accN, usa la routine con ghost-knots;
%    altrimenti usa la versione "standard" v1=vN=0.
%
% Esempi ------------------
%   quick_spline_demo(t, Q);                     % v1=vN=0
%   quick_spline_demo(t, Q, acc0, []);           % impone acc in t1
%   quick_spline_demo(t, Q, acc0, accN);         % acc in t1 e tN


    if nargin < 3 || isempty(method),  method = 'cubicC2';  end

    switch lower(method)
        case 'cubicc2'                          % le tue routine originali
            if nargin <= 3
                S = spline_cubic_tridiag(t,Q);
            else
                acc0 = varargin{1};
                accN = [];
                if nargin >= 4,  accN = varargin{2};  end
                S = spline_cubic_accbord(t,Q,acc0,accN);
            end

        case 'pchip'                            % nuovo ramo
            if nargin > 3
                warning('Acc0/AccN ignorate in P-chip');
            end
            S = spline_pchip(t,Q);

        otherwise
            error('Metodo sconosciuto');
    end

    % tabella per ogni giunto
    fprintf('\n=== VALORI AI NODI ===\n');
    [nJ, N] = size(Q);
    for j = 1:nJ
        T = table(t(:), Q(j,:).', S.v(j,:).', S.acc(j,:).', ...
                   'VariableNames',{'t','q','v','acc'});
        fprintf('\nGiunto %d\n', j), disp(T)
    end

    % grafico
    spline_utility_plot(t, Q, S);
end
