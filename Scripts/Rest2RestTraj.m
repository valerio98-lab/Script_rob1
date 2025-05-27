function [t,s,v,a] = Rest2RestTraj(L, A_max, V_max, varargin)
% MOTIONPROFILE  –  profili bang–coast–bang con V0 arbitrario
%
%   [t,s,v,a] = MotionProfile(L, A_max, V_max, 'Profile','State2Rest', ...
%                             'V0', -1.5, 'Scale',1.2)
%
% Name–Value:
%   'V0'      iniziale (default 0)   • segno reale
%   'Scale'   fattore temporale k>0  • default 1
%   'Profile' 'Rest2Rest' | 'State2Rest' (default 'Rest2Rest')
%
% Il segno di L stabilisce la direzione finale (+ccw, –cw).
% Non si ribaltano i grafici a posteriori: tutti i segni sono reali.
% -----------------------------------------------------------------------

% ---------------------------------------------------------------------
% TEORIA DI RIFERIMENTO
%   • Il profilo "bang-coast-bang" è composto da 3 fasi:
%       a) Accelerazione costante  +A_max   per Ts
%       b) Velocità costante       V_max    per Tc
%       c) Decelerazione costante  –A_max   per Ts
%   • Per profili "triangolari" Tc = 0 (ossia V_pk < V_max).
%   • Le relazioni fondamentali sono:
%          Ts = V_max / A_max                           (1)
%          L  = V_max·Tc + V_max² / A_max               (2)  (trapezoide)
%          L  = ½·A_max·Ts² + ½·A_max·Ts² = A_max·Ts²   (3)  (triangolo)
%   • Nel caso State-to-Rest la prima fase scompare:
%          L  = ½·V0² / A_max + V0·Tc                   (4)
%          Ts = V0 / A_max                              (5)
%   • Con la scalatura temporale k>1 (più lento) si ha:
%         a_scaled = a/k² ,    v_scaled = v/k ,    t_scaled = t·k
%
% ---------------------------------------------------------------------
%‒‒‒‒ parsing ------------------------------------------------------------
p = inputParser;
addParameter(p,'V0',0,@isnumeric);
addParameter(p,'Scale',1,@(x)x>0);
addParameter(p,'Profile','Rest2Rest', ...
             @(s)any(strcmpi(s,{'Rest2Rest','State2Rest'})));
parse(p,varargin{:});
V0   = p.Results.V0;
k    = p.Results.Scale;
mode = lower(p.Results.Profile);

sgn  = sign(L); if sgn==0, sgn = 1; end   % direzione target
Lmag = abs(L);                            % modulo spostamento
A    = A_max;  Vlim = V_max;

%‒‒‒ inizializza vettori globali (vuoti) -------------------------------
t = [];  a = [];  v = [];  s = [];

%======== nested function: concatena un segmento ========================
    function appendSeg(tt, aa, vv)
        % tt,aa,vv devono già essere vettori colonna (o riga stessa dim.)
        tt = tt(:);  aa = aa(:);  vv = vv(:);      % garantisco colonne
        if isempty(t)
            t = tt;  a = aa;  v = vv;
            s = cumtrapz(t,v);                     % integrazione da 0
        else
            tt = tt + t(end);                      % shift temporale
            t  = [t ; tt(2:end)];                  % niente punto doppio
            a  = [a ; aa(2:end)];
            v  = [v ; vv(2:end)];
            ds = cumtrapz(tt,vv);                  % integrazione locale
            s  = [s ; s(end) + ds(2:end)];
        end
    end
%========================================================================

% ①  frenata preliminare se V0 ha segno opposto a L
if strcmp(mode,'state2rest') && sign(V0) ~= 0 && sign(V0) ~= sgn
    Td   = abs(V0)/A;
    tt   = linspace(0,Td,200);
    aa   = -sign(V0)*A * ones(size(tt));
    vv   = V0 + aa.*tt;                           % v(t) in frenata
    appendSeg(tt, aa, vv);
    Lmag = Lmag - sgn*V0^2/(2*A);                 % tratto residuo
    V0   = 0;                                     % ora fermo
end

% ②  profilo bang–coast–bang sul tratto rimanente (partenza V0, arrivo 0)
V0 = min(max(V0,-Vlim),Vlim);                     % clamp a ±Vlim

Ts0 = (Vlim - V0)/A;                              % tempo → Vlim
if Lmag >= (Vlim^2 - V0^2)/(2*A)                  % TRAPEZOIDALE
    Ts  = Ts0;
    Tc  = (Lmag*A - (Vlim^2 - V0^2)/2)/(A*Vlim);
    Vpk = Vlim;
else                                              % TRIANGOLARE (Tc=0)
    Ts  = (-V0 + sqrt(V0^2 + 2*A*Lmag))/A;
    Tc  = 0;
    Vpk = V0 + A*Ts;
end

Ts = Ts*k;  Tc = Tc*k;                            % scalatura temporale

N   = 300;                                        % punti per segmento
% -- fase 1: accelerazione
t1  = linspace(0,Ts,N);
a1  =  sgn*A * ones(size(t1));
v1  =  V0 + sgn*A*t1;
appendSeg(t1,a1,v1);

% -- fase 2: cruise (se presente)
if Tc > 0
    t2 = linspace(0,Tc,N);
    a2 = zeros(size(t2));
    v2 =  Vpk*sgn * ones(size(t2));
    appendSeg(t2,a2,v2);
end

% -- fase 3: decelerazione
t3  = linspace(0,Ts,N);
a3  = -sgn*A * ones(size(t3));
v3  =  Vpk*sgn - sgn*A*t3;
appendSeg(t3,a3,v3);

% ③  plot se l’utente non chiede output
figure('Name',['MotionProfile – ' mode],'Units','normalized',...
       'Position',[.14 .23 .58 .58]);
subplot(3,1,1)
    plot(t,s,'b','LineWidth',1.4), grid on
    title('Posizione s(t)'), xlabel('t [s]')
subplot(3,1,2)
    plot(t,v,'r','LineWidth',1.4), grid on
    title('Velocità v(t)'),  xlabel('t [s]')
subplot(3,1,3)
    plot(t,a,'g','LineWidth',1.4), grid on
    title('Accelerazione a(t)'), xlabel('t [s]')
end




% ────────────────────────────────────────────────────────────────────────
% ❑   PROFILI NON NATIVAMENTE GESTITI   ❑
% ────────────────────────────────────────────────────────────────────────
%
%  REST → STATE   (parto da fermo, arrivo con Vf ≠ 0)
%  ─────────────
%  • Creo il profilo “speculare” STATE→REST invertendo tempo e segni:
%
%       % profilo di riferimento (state→rest con V0 = +Vf e stesso L)
%       [t0,s0,v0,a0] = MotionProfile(L, A, Vmax, ...
%                                     'Profile','State2Rest', ...
%                                     'V0',  Vf);
%
%       % ribalto il tempo (t = T - t0)  e cambio segno all’accelerazione
%       T   = t0(end);
%       t   = T - t0;
%       s   = flipud(s0);               % s(t) resta con stesso segno
%       v   = -flipud(v0);              % v(t) si ribalta di segno
%       a   =  flipud(a0);              % a(t) cambia verso con il ribaltamento
%
%  STATE → STATE  (V0 → Vf, entrambi ≠ 0)
%  ─────────────
%  • Spezzo il moto in due:
%        ① V0 → 0  (state→rest)  su un tratto L1
%        ② 0  → Vf (rest →state) sul tratto L2 = L - L1
%    L1 è la distanza percorsa per frenare da V0 a 0:
%        L1 = |V0|² / (2·A)
%
%       % fase-1: frenata
%       [t1,s1,v1,a1] = MotionProfile(sgn*L1, A, Vmax, ...
%                                     'Profile','State2Rest', ...
%                                     'V0', V0);
%
%       % fase-2: rest→state (ottenuta come sopra)
%       L2 = L - sgn*L1;                % sgn = sign(L)
%       [t2s,s2s,v2s,a2s] = MotionProfile(L2, A, Vmax, ...
%                                         'Profile','State2Rest', ...
%                                         'V0', abs(Vf));  % profilo speculare
%       T2 = t2s(end);
%       t2 = T2 - t2s;                  % ribalta tempo
%       s2 = s1(end) + flipud(s2s);     % posizioni concatenate
%       v2 = -flipud(v2s);
%       a2 =  flipud(a2s);
%
%       % concatena le due fasi (occhio a eliminare il campione doppio)
%       t = [t1 ; t1(end)+t2(2:end)];
%       s = [s1 ; s2(2:end)];
%       v = [v1 ; v2(2:end)];
%       a = [a1 ; a2(2:end)];
%
%  NB: nel collegare i due segmenti assicurati di:
%      • usare lo stesso segno di L (sgn) in entrambe le fasi
%      • rimuovere il primo campione della fase-2 per evitare sovrapposizioni
%
%  Con questi due “trucchi” puoi coprire tutti i casi (RR, SR, RS, SS)

% ────────────────────────────────────────────────────────────────────────
% Rest → State
% prendi il profilo State → Rest (già gestito) con V0 = +Vf.
% ribalti la sequenza temporale (t = T−t) → le curve vanno da 0 a +Vf.
% inverti i segni giusti (v cambia verso, a rimane coerente).

% State → State
% freni fino a fermarti (fase 1),
% poi usi il trucco precedente (fase 2) per ripartire fino a Vf.
% infine “cuci” i due segmenti togliendo il campione duplicato.
% Salvalo come commento: quando ti servirà basterà copiare-incollare quei pochi comandi nel tuo script di test.

