function S = naturalSpline(t, Q, varargin)
% naturalSpline  –  Cubic spline C² “natural” + grafici.

%% Per stampare i coefficienti in maniera leggibile: 
% [nJ, P, ~] = size(S.coeffs);       % P = N-1 = numero tratti
% for k = 1:P
%     fprintf('\n=== Tratto %d ===\n',k)
%     for j = 1:nJ
%         c = squeeze(coeffs(j,k,:));
%         fprintf('joint %d:  %.4f  %+ .4f τ  %+ .4f τ²  %+ .4f τ³\n', ...
%                 j, c(1), c(2), c(3), c(4));
%     end
% end
%%
%
%   S = naturalSpline(t, Q)                  % plot automatico
%   S = naturalSpline(t, Q, 'Samples',800)   % griglia più fitta
%   S = naturalSpline(t, Q, 'Plot',false)    % solo calcoli, niente plot
%
% INPUT
%   t : 1×N     nodi temporali strettamente crescenti
%   Q : nJ×N    valori (una riga per giunto)
%
% OPZIONI (Name–Value)
%   'Samples' : punti di campionamento (default 600)
%   'Plot'    : true/false (default true)
%
% OUTPUT (struct)
%   S.vNodes      nJ×N   velocità ai nodi
%   S.coeffs      nJ×(N-1)×4   [a0 a1 a2 a3] per tratto
%   S.tGrid       1×M    griglia
%   S.q / .v / .a nJ×M   pos-vel-acc campionate
%
% ────────────────────────────────────────────────────────────────
% Usa la funzione cubicCoeffs(*) che possiedi già per costruire
% i coefficienti di ogni tratto; qui la richiamo solo per i 4
% parametri [a0 a1 a2 a3] sfruttando le velocità ai nodi
% calcolate con il sistema tridiagonale “natural”.
%
% (*) se la tua cubicCoefficients si chiama diverso, cambia 
%     la riga marcata con >>> in fondo al file.

% ------------------------------------------------– parsing
p  = inputParser;
addParameter(p,'Samples',600,@(x)isnumeric(x)&&isscalar(x)&&x>10);
addParameter(p,'Plot',   true,@islogical);
parse(p,varargin{:});
M       = p.Results.Samples;
doPlot  = p.Results.Plot;

t = t(:).';                          % assicura vettore riga
[nJ,N] = size(Q);
if numel(t)~=N
    error('Dimensioni incompatibili fra t e Q');
end

% ------------------------------------------------– velocità ai nodi (bordo NATURAL)
vNodes = zeros(nJ,N);
for j = 1:nJ
    vNodes(j,:) = cubicSplineVelNatural(t, Q(j,:));
end


% ---------- coefficienti cubici in tempo fisico -------------------
h   = diff(t);
coeffs = zeros(nJ,N-1,4);

for k = 1:N-1
    h_k = h(k);
    qi  = Q(:,k);      qf = Q(:,k+1);
    vi  = vNodes(:,k); vf = vNodes(:,k+1);

    a0 = qi;
    a1 = vi;
    a2 = ( 3*(qf-qi)/h_k - 2*vi - vf ) / h_k;
    a3 = ( vi + vf - 2*(qf-qi)/h_k )   / h_k^2;

    coeffs(:,k,:) = [a0 a1 a2 a3];
end


% ------------------------------------------------– campionamento su griglia
tGrid = linspace(t(1), t(end), M);
[qGrid, vGrid, aGrid] = evalSplineGrid(tGrid, t, coeffs);

% ------------------------------------------------– output struct
S.vNodes  = vNodes;
S.coeffs  = coeffs;
S.tGrid   = tGrid;
S.q       = qGrid;
S.v       = vGrid;
S.a       = aGrid;

format short g                 % meno zeri inutili
for k = 1:size(coeffs,2)       % ogni tratto
    fprintf('\nTratto %d (h = %.3f s):\n',k, diff(t(k:k+1)))
    for j = 1:size(coeffs,1)
        a = squeeze(coeffs(j,k,:));
        fprintf('  joint %d:  %.4f  %+ .4f τ  %+ .4f τ²  %+ .4f τ³\n',...
                j,a(1),a(2),a(3),a(4));
    end
end


% ------------------------------------------------– plot
if doPlot
    col = lines(nJ);    lg = arrayfun(@(j) sprintf('joint %d',j),1:nJ,'uni',0);
    figure('Name','Natural cubic spline'), tiledlayout(3,1)
    nexttile, hold on, grid on, title('q(t)')
    for j=1:nJ, plot(tGrid,qGrid(j,:), 'Color',col(j,:)), ...
                   plot(t,Q(j,:),'o','Color',col(j,:),'MarkerFace','w'); end
    legend(lg), ylabel('rad')

    nexttile, hold on, grid on, title('q̇(t)')
    for j=1:nJ, plot(tGrid,vGrid(j,:), 'Color',col(j,:)); end
    ylabel('rad/s')

    nexttile, hold on, grid on, title('q̈(t)')
    for j=1:nJ, plot(tGrid,aGrid(j,:), 'Color',col(j,:)); end
    ylabel('rad/s²'), xlabel('time  [s]')
end
end  % ===== end naturalSpline =======================================


% ===== helper – velocità ai nodi con condizione NATURAL ============
function v = cubicSplineVelNatural(t,q)
t = t(:).';  q = q(:).';
N = numel(t);   h = diff(t);   d = diff(q)./h;

% matrice tridiagonale per nodi interni
A   = spdiags([[h(2:end-1),0]'  2*(h(1:end-1)+h(2:end))'  [0,h(2:end-1)]'], ...
              [-1 0 1], N-2, N-2);
rhs = 3*(h(1:end-1).*d(2:end) + h(2:end).*d(1:end-1)).';

v = zeros(1,N);
if N>2
    v(2:N-1) = (A\rhs).';
end
v(1)   = (3*d(1)   - v(2))   / 2;
v(end) = (3*d(end) - v(end-1))/ 2;
end

% ===== helper – valuta spline su griglia ===========================
function [qG,vG,aG] = evalSplineGrid(tt,t, C)
nJ = size(C,1);     N = numel(t);
qG = zeros(nJ,numel(tt));    vG = qG;    aG = qG;

for k=1:N-1
    idx = tt>=t(k) & tt<=t(k+1);
    tau = tt(idx) - t(k);

    Tpos = [ones(size(tau)); tau; tau.^2; tau.^3];
    for j=1:nJ
        c = squeeze(C(j,k,:));
        qG(j,idx) = c.'*Tpos;
        vG(j,idx) = c(2) + 2*c(3)*tau + 3*c(4)*tau.^2;
        aG(j,idx) = 2*c(3) + 6*c(4)*tau;
    end
end
end
