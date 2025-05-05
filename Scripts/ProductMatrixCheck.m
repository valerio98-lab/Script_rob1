function matrix = ProductMatrixCheck(varargin)
% Product of an arbitrary number of rotation matrices.
% Usage examples:
%   Rif = ProductMatrixCheck(R1, R2);                       % solo matrici
%   Rif = ProductMatrixCheck(R1, R2, true, r, theta);       % con controllo di Rodriguez

    nArg = nargin;
    assert(nArg >= 1, 'Serve almeno una matrice di input.');

    % --- Riconosci se l’utente ha chiesto il controllo Rodriguez -----------
    checkRodriguez = false;
    if  islogical(varargin{end})                   % … ultimo argomento è logical
        checkRodriguez = varargin{end};
        varargin(end) = [];                        % rimuovilo
    end

    if  checkRodriguez
        assert(numel(varargin) >= 3, ...
               'Per checkRodriguez=true servono anche "r" e "theta".');
        theta = varargin{end};   varargin(end) = [];
        r     = varargin{end};   varargin(end) = [];
    end
    % -----------------------------------------------------------------------

    % --- Prodotto delle matrici --------------------------------------------
    prodMatrix = eye(3);
    for k = 1:numel(varargin)
        prodMatrix = prodMatrix * varargin{k};
    end

    % -----------------------------------------------------------------------

    % --- Eventuale verifica -------------------------------------------------
    if checkRodriguez
        check_rotation_matrix(prodMatrix, true, theta, r);
    else
        check_rotation_matrix(prodMatrix);
    end

    matrix = prodMatrix;
end
