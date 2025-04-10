function R_finale = calcRotationMatrix(alpha, beta, gamma, order)

    % Verifica che l'ordine abbia esattamente 3 caratteri
    if length(order) ~= 3
        error('La stringa dell''ordine deve contenere esattamente 3 caratteri (X, Y o Z).');
    end
    
    % Definizione delle funzioni anonime per le rotazioni attorno agli assi
    RotX = @(theta) [1, 0, 0; 
                     0, cos(theta), -sin(theta); 
                     0, sin(theta), cos(theta)];
                 
    RotY = @(theta) [cos(theta), 0, sin(theta); 
                     0, 1, 0; 
                     -sin(theta), 0, cos(theta)];
                 
    RotZ = @(theta) [cos(theta), -sin(theta), 0; 
                     sin(theta), cos(theta), 0; 
                     0, 0, 1];
    
    % Array degli angoli per una facile indicizzazione
    angles = [alpha, beta, gamma];
    
    % Inizializza R_finale con la matrice identità
    R_finale = eye(3);
    
    % Applica le rotazioni in base all'ordine specificato
    for i = 1:3
        switch upper(order(i))
            case 'X'
                R_finale = R_finale * RotX(angles(i));
            case 'Y'
                R_finale = R_finale * RotY(angles(i));
            case 'Z'
                R_finale = R_finale * RotZ(angles(i));
            otherwise
                error('Asse non valido: %s. Usa solo X, Y o Z.', order(i));
        end
    end

end