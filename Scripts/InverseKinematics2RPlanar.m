function [allSolutions, numSol] = InverseKinematics2RPlanar(l1,l2, px, py, q1_min, q1_max, q2_min, q2_max)
% IK2RPlanar Restituisce TUTTE le soluzioni inverse (q1,q2) di un robot 2R planare 
% per la posizione cartesiana (px,py).
%
% Parametri:
%   l1, l2: lunghezze dei link
%   px, py: coordinate finali del punto da raggiungere
%
% Output:
%   allSolutions: matrice di dimensioni (numSol x 2) contenente tutte le 
%                 soluzioni [q1, q2] trovate. 
%                 (se numSol=2, ci saranno due righe)
%   numSol: numero di soluzioni trovate (0,1,2).

    % Calcoliamo cos(q2)
    if nargin < 8 || isempty(q2_max)
        q2_max = +inf;
    end
    if nargin < 7 || isempty(q2_min)
        q2_min = -inf;
    end
    if nargin < 6 || isempty(q1_max)
        q1_max = +inf;
    end
    if nargin < 5 || isempty(q1_min)
        q1_min = -inf;
    end

    c2 = (px^2 + py^2 - (l1^2 + l2^2)) / (2*l1*l2);
    
    % Se c2 è fuori dall'intervallo [-1, 1], nessuna soluzione
    if abs(c2) > 1
        allSolutions = [];     % nessuna soluzione
        numSol = 0;
        return
    end
    
    % Se c2 = ±1 (entro una piccola tolleranza numerica), abbiamo 1 soluzione
    epsilon = 1e-12;
    if abs(abs(c2) - 1) < epsilon
        % In questo caso, s2 ~ 0 => q2 è 0 o pi (se c2=-1).
        % s2 = 0 => atan2(0, ±1) = 0 o pi
        s2 = 0;
        q2 = atan2(s2, c2);
        
        % Calcoliamo q1
        % q1 = atan2(py, px) - atan2(l2 * s2, l1 + l2 * c2)
        q1 = atan2(py, px) - atan2(l2*s2, l1 + l2*c2);
        
        allSolutions = [q1, q2];
        numSol = 1;
        
    else
        % Altrimenti, c2 è strettamente fra -1 e 1, quindi ci sono 2 soluzioni
        s2_pos =  sqrt(1 - c2^2);   % gomito in alto
        s2_neg = -sqrt(1 - c2^2);   % gomito in basso
        
        q2_pos = atan2(s2_pos, c2);
        q2_neg = atan2(s2_neg, c2);
        
        % Calcoliamo i corrispondenti q1
        q1_pos = atan2(py, px) - atan2(l2*s2_pos, l1 + l2*c2);
        q1_neg = atan2(py, px) - atan2(l2*s2_neg, l1 + l2*c2);
        
        solutions = [q1_pos, q2_pos;
                        q1_neg, q2_neg];
        numSol = 2;
        
        allSolutions = [];
        for i = 1:numSol
            q1 = solutions(i,1);
            q2 = solutions(i,2);
            if (q1 >= q1_min && q1 <= q1_max) && (q2 >= q2_min && q2 <= q2_max)
                allSolutions = [allSolutions; q1, q2];
            else 
                numSol = numSol - 1; 
            end
        end


    end
    
end
