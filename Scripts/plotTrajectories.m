function plotTrajectories(a0,a1,a2,a3,T)
    % a0..a3 sono vettori n x 1 dei coefficienti per ciascun joint
    % T è il tempo totale di percorrenza
    n = numel(a0);
    
    % costruiamo un campionamento in "tempo reale"
    N = 200;
    t = linspace(0, T, N);      % tempo [0 .. T]
    s = t / T;                  % parametro normalizzato [0 .. 1]
    
    % preallocazione
    Q   = zeros(n, N);
    Qd  = zeros(n, N);
    
    % per ogni joint i calcoliamo posizione e velocità
    for i = 1:n
        Q(i,:)  = a0(i) + a1(i)*s + a2(i)*s.^2 + a3(i)*s.^3;
        Qd(i,:) =        a1(i) ...
                  + 2*a2(i)*s   ...
                  + 3*a3(i)*s.^2;      % derivata rispetto a s
        Qd(i,:) = Qd(i,:) / T;        % d/dt = (d/ds)*(ds/dt) = (…)/T
    end
    
        % plot posizioni
    % Q(i,:)  = posizione del joint i in funzione di t
    % Qd(i,:) = velocità    del joint i in funzione di t
    % t = vettore dei tempi 1×N
    
    % Colori “automatici” per distinguerli
    colors = lines(n);
    
    % ---- Sovrappongo le posizioni ----
    figure; hold on;
    for i = 1:n
        plot(t, Q(i,:), 'LineWidth',1.5, 'Color', colors(i,:));
    end
    xlabel('Time (s)');
    ylabel('Position (rad)');
    title('Overlay Joint Positions');
    legend(arrayfun(@(i)sprintf('q_%d',i),1:n,'UniformOutput',false), ...
           'Location','best');
    grid on;
    hold off;
    
    % ---- Sovrappongo le velocità ----
    figure; hold on;
    for i = 1:n
        plot(t, Qd(i,:), 'LineWidth',1.5, 'Color', colors(i,:));
    end
    xlabel('Time (s)');
    ylabel('Velocity (rad/s)');
    title('Overlay Joint Velocities');
    legend(arrayfun(@(i)sprintf('q_%d dot',i),1:n,'UniformOutput',false), ...
           'Location','best');
    grid on;
    hold off;

end
