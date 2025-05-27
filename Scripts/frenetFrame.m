function FF = frenetFrame(p, s)
% FRENETFRAME Calcola simbolicamente il Frenet–Serret frame
%
%   FF = frenetFrame(p, s) prende in input
%     p : 3×1 vettore simbolico funzione di s, ad es. [x(s); y(s); z(s)]
%     s : simbolico, variabile di parametrizzazione
%   e restituisce una struct FF con i campi:
%     FF.T      : versore tangente   T(s)
%     FF.N      : versore normale   N(s)
%     FF.B      : versore binormale  B(s)
%     FF.kappa  : curvatura κ(s)
%     FF.tau    : torsione  τ(s)
%
%   Usa le definizioni:
%     T = p′/||p′||
%     B = (p′ × p″) / ||p′ × p″||
%     N = B × T
%     κ = ||p′ × p″|| / ||p′||^3
%     τ = [p′ · (p″ × p‴)] / ||p′ × p″||^2

    % derivate prime, seconde e terze
    p1 = diff(p, s);    % p′
    p2 = diff(p1, s);   % p″
    p3 = diff(p2, s);   % p‴

    % tangente unitario
    T = simplify((p1/norm(p1)),"IgnoreAnalyticConstraints",true);

    % normale
    N = (cross(p1, cross(p2,p1)))/(norm(p1)*norm(cross(p2,p1)));
    N = simplify(N, 'IgnoreAnalyticConstraints',true); 

    % normale
    B = simplify( cross(T, N), 'IgnoreAnalyticConstraints',true);

    % curvatura
    kappa = simplify( norm(cross(p1,p2)) / norm(p1)^3, "IgnoreAnalyticConstraints",true );

    % torsione
    tau = simplify(dot(p1,cross(p2,p3)) / norm(cross(p1,p2))^2, "IgnoreAnalyticConstraints",true );

    % componi la struct di output
    FF = struct( ...
      'T',      T, ...
      'N',      N, ...
      'B',      B, ...
      'kappa',  kappa, ...
      'tau',    tau ...
    );

    % (facoltativo) stampa a video
    % disp('T(s) = '); pretty(T)
    % disp('N(s) = '); pretty(N)
    % disp('B(s) = '); pretty(B)
    % disp('κ(s) = '); pretty(kappa)
    % disp('τ(s) = '); pretty(tau)
end
