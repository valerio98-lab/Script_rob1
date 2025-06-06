function output = BinaryGreyDecimal(value, fromType, toType)
% convertCode  Converte tra Gray, binario e decimale.
%
%   OUTPUT = convertCode(VALUE, FROMTYPE, TOTYPE) accetta:
%     - VALUE: stringa di '0' e '1' se FROMTYPE è 'gray' o 'binary', 
%              oppure un numero intero non negativo se FROMTYPE è 'decimal'.
%     - FROMTYPE e TOTYPE: ognuno può essere 'gray', 'binary' o 'decimal'.
%
%   Esempi:
%     convertCode('1101', 'gray', 'binary')
%     convertCode('1011', 'binary', 'decimal')
%     convertCode(13, 'decimal', 'binary')
%     convertCode(13, 'decimal', 'gray')
%     convertCode('1110', 'binary', 'gray')
%
%   Se FROMTYPE == TOTYPE, restituisce VALUE così com’è.

    % Normalizzo i parametri in minuscolo
    fromType = lower(fromType);
    toType   = lower(toType);

    % Se il formato di partenza è uguale a quello di arrivo, restituisco direttamente
    if strcmp(fromType, toType)
        output = value;
        return;
    end

    % Passaggio intermedio: portare tutto in binario (stringa), poi da binario 
    % passare a destinazione (decimal o gray).
    switch fromType
        case 'gray'
            % VALUE è stringa di Gray: converto in binario (stringa)
            if ~ischar(value) && ~isstring(value)
                error('Per fromType = ''gray'' il valore deve essere una stringa di bit.');
            end
            binStr = gray2binary(char(value));
        case 'binary'
            % VALUE è stringa di bit: la tengo come binStr
            if ~ischar(value) && ~isstring(value)
                error('Per fromType = ''binary'' il valore deve essere una stringa di bit.');
            end
            binStr = char(value);
        case 'decimal'
            % VALUE è un numero intero non negativo
            if ~isscalar(value) || value < 0 || value ~= floor(value)
                error('Per fromType = ''decimal'' il valore deve essere un intero non negativo.');
            end
            % Converto in stringa binaria senza prefissi, esattamente sufficienti bit
            binStr = dec2bin(value);
        otherwise
            error('FROMTYPE deve essere ''gray'', ''binary'' o ''decimal''.');
    end

    % Ora binStr contiene sempre la rappresentazione in base 2 (stringa).
    % Passo quindi a TOTYPE:
    switch toType
        case 'binary'
            output = binStr;

        case 'decimal'
            % Converte la stringa binaria in intero
            output = bin2dec(binStr);

        case 'gray'
            % Converte la stringa binaria in code di Gray
            output = binary2gray(binStr);

        otherwise
            error('TOTYPE deve essere ''gray'', ''binary'' o ''decimal''.');
    end
end

% -------------------------------------------------------------------------
function binStr = gray2binary(grayStr)
% gray2binary  Converte una stringa di bit in Gray in stringa di bit binaria.
%
%   binStr = gray2binary(grayStr) assume grayStr come char array di '0' e '1'.
%   Restituisce binStr, stringa binaria equivalente (stessa lunghezza).

    n = length(grayStr);
    binVec = zeros(1, n);  % vettore di bit per il risultato binario
    % Primo bit: binario = stesso bit di Gray
    binVec(1) = str2double(grayStr(1));

    % Ogni bit successivo: bin(i) = bin(i-1) XOR gray(i)
    for i = 2:n
        g = str2double(grayStr(i));
        binVec(i) = xor(binVec(i-1), g);
    end

    % Ricostruisco la stringa da binVec
    binStr = char(binVec + '0');
end

% -------------------------------------------------------------------------
function grayStr = binary2gray(binStr)
% binary2gray  Converte una stringa di bit binaria in stringa di bit in Gray.
%
%   grayStr = binary2gray(binStr) assume binStr come char array di '0' e '1'.
%   Restituisce grayStr, stringa di Gray equivalente (stessa lunghezza).

    n = length(binStr);
    grayVec = zeros(1, n);
    binVec  = binStr - '0';  % converto in vettore numerico 0/1

    % Primo bit di Gray = primo bit di binario
    grayVec(1) = binVec(1);

    % Ogni bit successivo: gray(i) = bin(i-1) XOR bin(i)
    for i = 2:n
        grayVec(i) = xor(binVec(i-1), binVec(i));
    end

    % Ricostruisco la stringa da grayVec
    grayStr = char(grayVec + '0');
end
