function spline_utility_plot(t, Q, S, nGrid)
% plot_spline_multi  –  visualizza la traiettoria generata da
%                       spline_cubic_tridiag  **o**
%                       spline_cubic_accbord
%
%   t  : 1×N      (nodi reali, gli stessi passati allo solver)
%   Q  : nJ×N     (posizioni ai nodi reali)
%   S  : struct   (campo .C con i coefficienti cubici)
%   nGrid (opt) : numero di punti del grafico – default 600
%
%   NB: la funzione riconosce automaticamente se i coefficienti in
%       S.C includono o meno i tratti fantasma e prende quelli giusti.

    if nargin < 4, nGrid = 600; end
    [nJ, N] = size(Q);

    % 1) ricostruisci le posizioni su griglia fine
    tt = linspace(t(1), t(end), nGrid);
    qq = zeros(nJ, numel(tt));

    % C può contenere tratti fantasma → allinea con t (N-1 reali)
    C = S.C(:, 1:N-1, :);

    for k = 1:N-1
        idx = tt >= t(k) & tt <= t(k+1);
        tau = tt(idx) - t(k);
        for j = 1:nJ
            coeff = squeeze(C(j,k,:));          % [a0 a1 a2 a3]
            qq(j,idx) = polyval(flipud(coeff), tau);
        end
    end

    % 2) disegno
    col = lines(nJ);
    h = gobjects(1, nJ); 

    figure('Name','Cubic C^2 spline'), clf, hold on, grid on
    for j = 1:nJ
        h(j) = plot(tt, qq(j,:), 'Color', col(j,:), 'LineWidth',1.5, ...
            'DisplayName', sprintf('joint %d', j)); 
        plot(t,  Q(j,:),  'o', 'Color', col(j,:), ...
             'MarkerFaceColor','w', 'HandleVisibility', 'off');
    end
    legend(h, 'Location','best'); 
end
