function [dphi, T] = phiDotRotational(omega, phi, seq)
% phiDotGeneral  Calcola dphi per qualunque sequenza di rotazioni fisse
%
% INPUT:
%   R     : matrice di rotazione 3×3 (numeric)
%   omega : vettore angolar-velocità [ωx;ωy;ωz]
%   phi   : vettore angoli [φ1;φ2;φ3] in radianti
%   seq   : stringa di 3 caratteri tra 'x','y','z', es. 'xzy'
%
% OUTPUT:
%   dphi  : vettore [φ̇1; φ̇2; φ̇3]
%   T     : matrice 3×3 tale che ω = T * φ̇

    % 1) Definiamo simbolicamente gli angoli e costruiamo R_sym
    syms p1 p2 p3 real
    % map lettere → rotazioni di assi fissi
    Rx = @(a)[1 0         0;
              0 cos(a) -sin(a);
              0 sin(a)  cos(a)];
    Ry = @(a)[ cos(a) 0 sin(a);
               0      1 0;
              -sin(a) 0 cos(a)];
    Rz = @(a)[cos(a) -sin(a) 0;
              sin(a)  cos(a) 0;
              0       0      1];
    % costruiamo la sequenza (extrinsic fixed axes)
    R_sym = eye(3);
    for i = 1:3
        switch lower(seq(i))
          case 'x', R_sym = R_sym * Rx([p1 p2 p3](i));
          case 'y', R_sym = R_sym * Ry([p1 p2 p3](i));
          case 'z', R_sym = R_sym * Rz([p1 p2 p3](i));
          otherwise, error('Sequenza invalida');
        end
    end

    % 2) calcolo le derivate parziali e riempio T_sym
    T_sym = sym(zeros(3));
    R_T_sym = transpose(R_sym);
    for i = 1:3
        dR_dpi = diff(R_sym, [p1 p2 p3](i));        % ∂R/∂φᵢ
        S_i    = simplify(dR_dpi * R_T_sym);         % S_i = ∂R/∂φᵢ * Rᵀ
        % vee-operator: estraggo [S(3,2); S(1,3); S(2,1)]
        T_sym(:,i) = [ S_i(3,2); S_i(1,3); S_i(2,1) ];
    end

    % 3) sostituisco i valori numerici di φ nella T_sym
    T = double( subs(T_sym, [p1 p2 p3], phi(:).') );

    % 4) risolvo T * dphi = omega
    dphi = T \ omega;
end
