function [v_dir, sigma, cosines] = AnalyzeJointContribution(J, mode)
% ANALYZEJOINTCONTRIBUTION 
%   [v_dir, sigma, cosines] = analyzeJointContribution(J, mode)
%
% Calcola tramite SVD la direzione di giunto che che massimizza/minimizza la velocità in uscita:
%   - massimizza ||J*v_dir|| (mode = 'max')
%   - minimizza ||J*v_dir|| (mode = 'min')
%
% Input:
%   J    : m×n Jacobian numerico
%   mode : 'max' (default) o 'min'
%
% Output:
%   v_dir   : n×1 vettore unitario in joint‐space
%   sigma   : singolo valore singolare corrispondente
%   cosines : n×1 array dei coseni tra v_dir ed e_i (joint axis)
%
    if nargin<2, mode = 'max'; end

    % calcolo SVD di J
    [~, S, V] = svd(J);

    switch lower(mode)
      case 'max'
        sigma = S(1,1);
        v_dir = V(:,1);
      case 'min'
        sigma = S(end,end);
        v_dir = V(:,end);
      otherwise
        error('mode deve essere ''max'' o ''min''');
    end

    % V è ortonormale => v_dir ha norma 1, i suoi elementi sono i coseni
    cosines = v_dir;

    % visualizzazione sintetica
    fprintf('--- %s gain ---\n', upper(mode));
    fprintf('σ = %g\n', sigma);
    for i = 1:length(v_dir)
      fprintf(' q%-2d: cosθ = % .4f\n', i, cosines(i));
    end
end
