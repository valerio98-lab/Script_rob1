function [q_out, guesses, fk_history, cartesian_errors, alpha_hist] = GradientMethod(q_sym, desired_point, direct_kyn, q0, varargin)
% GRADIENTMETHOD   Inverse kinematics via pure (optionally damped) gradient descent
%

%% ESEMPI
% f_q = [0.5*cos(q1) + 0.5*cos(q1+q2)*cos(q3);
%     0.5*sin(q1)+0.5*sin(q1+q2)*cos(q3);
%     0.5+0.5*sin(q3)];


%   [q_out, guesses, fk_hist, errs, alpha_hist] = GradientMethod(q_sym, desired_point,
%                                     direct_kyn, q0, ...
%                                     'Verbose',true, ...
%                                     'Alpha',0.05, ...
%                                     'LineSearch',true, ...
%                                     'SkewGain',Ks, ... )
% 
% 
% [q_out, guesses, fk_history, cartesian_errors] = GradientMethod([q1 q2 q3], [0.3;-0.3;0.7], f_q, [-pi/4 pi/4 pi/4], ...
%     'CartesianTolerance', 10e-3,'Verbose',true, 'StoreTraj', true); 
%
%%
% INPUTS (positional)
%   q_sym          : 1×n symbolic row vector of joints (e.g. [q1 q2 ...])
%   desired_point  : m×1 numeric vector (target pose in task-space)
%   direct_kyn     : m×1 symbolic column  f(q_sym)   (forward kinematics)
%   q0             : 1×n numeric row (initial guess)
%
% OPTIONAL NAME–VALUE PAIRS
%   'Verbose'      (false)   if true, print iteration log
%   'MaxIterations'(1000)    maximum iterations
%   'CartesianTolerance' (1e-8) stop when ||e|| below threshold
%   'Alpha'        (0.05)    initial step size
%   'LineSearch'   (false)   enable backtracking line search on alpha
%   'Beta'         (0.5)     reduction factor for line search (0<β<1)
%   'AlphaMin'     (1e-6)    minimum alpha allowed in line search
%   'SkewGain'     ([])      m×m skew-symmetric matrix K  (uses (I+K)·e )
%   'StoreTraj'    (false)   if true, saves q at each iteration
%
% OUTPUTS
%   q_out          : n×1 solution
%   guesses        : iter×n matrix of joint guesses          (empty if StoreTraj=false)
%   fk_history     : m×iter matrix of f(q) at each step
%   cartesian_errors: 1×iter vector of ||e|| history
%   alpha_hist     : 1×iter vector of alpha actually used
%
% NOTES ---------------------------------------------------------------
% – Uses analytic Jacobian derived once from symbolic input.
% – Pure gradient step:  dq = alpha * J' * (I+K) * e
%   (if SkewGain is empty, K=0).
% – Optional simple backtracking line‐search guarantees monotone error
%   decrease.
% --------------------------------------------------------------------

% --- parse inputs ----------------------------------------------------
import matlab.internal.inputParser.*

p = inputParser;
addRequired (p,'q_sym',         @(x) isa(x,'sym') && isrow(x));
addRequired (p,'desired_point', @(x) isnumeric(x) && isvector(x));
addRequired (p,'direct_kyn',    @(x) isa(x,'sym'));
addRequired (p,'q0',            @(x) isnumeric(x) && isvector(x));
addParameter(p,'Verbose', false,       @(x) islogical(x)&&isscalar(x));
addParameter(p,'MaxIterations', 1000,  @(x) isnumeric(x)&&isscalar(x));
addParameter(p,'CartesianTolerance',1e-8,@(x) isnumeric(x)&&isscalar(x));
addParameter(p,'Alpha', 0.05,          @(x) isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'LineSearch', false,    @(x) islogical(x)&&isscalar(x));
addParameter(p,'Beta', 0.5,            @(x) isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'AlphaMin', 1e-6,       @(x) isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'SkewGain', [],         @(x) isempty(x) || isnumeric(x));
addParameter(p,'StoreTraj', false,     @(x) islogical(x)&&isscalar(x));
parse(p,q_sym,desired_point,direct_kyn,q0,varargin{:});
opts = p.Results;

% --- build numeric handles ------------------------------------------
q_sym = q_sym(:).';
direct_kyn = reshape(direct_kyn,[],1);
f_num = matlabFunction(direct_kyn,'Vars',{q_sym});
J_sym = jacobian(direct_kyn,q_sym);
J_num = matlabFunction(J_sym,'Vars',{q_sym});

nd = numel(desired_point);
nq = numel(q_sym);
qi = opts.q0(:);
M  = opts.MaxIterations;

% sanity check on SkewGain -------------------------------------------
if isempty(opts.SkewGain)
    useK = false;
else
    K = opts.SkewGain;
    if size(K,1)~=nd || size(K,2)~=nd
        error('SkewGain must be a %d×%d matrix.',nd,nd);
    end
    if max(abs(K+K.'))>1e-10
        warning('SkewGain is not perfectly skew-symmetric (K^T ≠ -K).');
    end
    useK = true;
end

% allocate logs -------------------------------------------------------
maxLog = M;
if opts.StoreTraj, guesses = zeros(maxLog,nq); else, guesses = []; end
fk_history       = zeros(nd,maxLog);
cartesian_errors = zeros(1,maxLog);
alpha_hist       = zeros(1,maxLog);

alpha = opts.Alpha;

% --- main loop -------------------------------------------------------
for k = 1:M
    p_k = f_num(qi.');          % current EE pose
    e_k = desired_point(:) - p_k;
    if useK
        e_mod = (eye(nd)+K)*e_k;
    else
        e_mod = e_k;
    end

    Jk  = J_num(qi.');          % Jacobian
    grad = Jk.' * e_mod;        % gradient direction

    % ---- optional backtracking line search -------------------------
    if opts.LineSearch
        alpha_ls = alpha;
        while alpha_ls >= opts.AlphaMin
            q_try = qi + alpha_ls * grad;
            err_try = desired_point(:) - f_num(q_try.');
            if norm(err_try) < norm(e_k) - 1e-12
                break;
            end
            alpha_ls = opts.Beta * alpha_ls;
        end
        if alpha_ls < opts.AlphaMin
            warning('Line search failed to improve. Using minimum alpha.');
            alpha_ls = opts.AlphaMin;
        end
        alpha_used = alpha_ls;
    else
        alpha_used = alpha;
    end

    % ---- update configuration -------------------------------------
    qi = qi + alpha_used * grad;

    % ---- save logs -------------------------------------------------
    fk_history(:,k)    = p_k;
    cartesian_errors(k)= norm(e_k);
    alpha_hist(k)      = alpha_used;
    if opts.StoreTraj, guesses(k,:) = qi.'; end

    % ---- stopping criteria ----------------------------------------
    if cartesian_errors(k) < opts.CartesianTolerance
        break;
    end
end

it = k;  % actual iterations executed
fk_history       = fk_history(:,1:it);
cartesian_errors = cartesian_errors(1:it);
alpha_hist       = alpha_hist(1:it);
if opts.StoreTraj
    guesses = guesses(1:it,:);
end

q_out = qi;

% --- verbose log -----------------------------------------------------
if opts.Verbose
    fprintf('\n--- GradientMethod history (iter %d) ---\n',it);
    if opts.StoreTraj
        for i = 1:it
            fprintf('Iter %3d: q=[%s],  err=%.3e,  α=%.3e\n', ...
                i, num2str(guesses(i,:),'%.6f '), ...
                cartesian_errors(i), alpha_hist(i) );
        end
    else
        for i = 1:it
            fprintf('Iter %3d: err=%.3e,  α=%.3e\n', ...
                i, cartesian_errors(i), alpha_hist(i) );
        end
    end
    fprintf('------------------------------------------\n');
end
end
