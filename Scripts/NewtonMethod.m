function [q_out, guesses, fk_history, cartesian_errors] = NewtonMethod(q_sym, desired_point, direct_kyn, q0, varargin)
% NEWTONMETHOD   Inverse kinematics via (damped) Newton iteration
%
% [q_out, guesses, fk_hist, errs] = NewtonMethod(q_sym, desired_point, direct_kyn, q0, ...
%                                   'Verbose, true, ...
%                                   'Method', 'dls')
%
% INPUTS (positional)
%   q_sym           : symbolic row vector of joints, e.g. [q1 q2 q3]
%   desired_point   : numeric column [px; py; pz] (o 2×1)
%   direct kyn      : symbolic column f(q_sym)
%   q0              : initial guess, numeric row [q10, q20, q30]
%
% OPTIONAL NAME–VALUE PAIRS
%    verbose                  (default false) %if true print all the values
%   'MaxIterations'           (default 1000)
%   'CartesianTolerance'      (default 1e-8)
%   'MinJointIncrement'       (default 1e-6)
%   'MaxClosenessSingularity' (default 1e-6)
%   'Method'                  {'auto','dls'} (default 'auto')
%   'DLSLambda'               (default 1e-2)    % solo se Method='dls'
%
% OUTPUTS
%   q_out           : solution (column)
%   guesses         : [iter×ndof] joint history
%   fk_history      : [nd×iter] f(q) history
%   cartesian_errors: [1×iter] error history

  %— parse inputs ---------------------------------------------------------
  p = inputParser;
  addRequired (p, 'q_sym',          @(x) isa(x,'sym') && isvector(x));
  addRequired (p, 'desired_point',  @(x) isnumeric(x) && isvector(x));
  addRequired (p, 'direct_kyn',          @(x) isa(x,'sym'));
  addRequired (p, 'q0',             @(x) isnumeric(x) && isvector(x));
  addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
  addParameter(p, 'MaxIterations',           1000,    @(x) isnumeric(x)&&isscalar(x));
  addParameter(p, 'CartesianTolerance',      1e-8,    @(x) isnumeric(x)&&isscalar(x));
  addParameter(p, 'MinJointIncrement',       1e-6,    @(x) isnumeric(x)&&isscalar(x));
  addParameter(p, 'MaxClosenessSingularity', 1e-6,    @(x) isnumeric(x)&&isscalar(x));
  addParameter(p, 'Method',            'auto',    @(x) ismember(lower(x),{'auto','dls'}));
  addParameter(p, 'DLSLambda',          1e-2,    @(x) isnumeric(x)&&isscalar(x));
  parse(p, q_sym, desired_point, direct_kyn, q0, varargin{:});
  opts = p.Results;

  %— prepare numeric functions -------------------------------------------
  q_sym = q_sym(:).';                    
  direct_kyn = reshape(direct_kyn, [], 1);         
  f_num = matlabFunction(direct_kyn, 'Vars',{q_sym});
  J_sym = jacobian(direct_kyn, q_sym);        
  J_num = matlabFunction(J_sym, 'Vars',{q_sym});

  nd = numel(desired_point);
  nq = numel(q_sym);
  qi = opts.q0(:);
  M  = opts.MaxIterations;

  guesses          = zeros(M, nq);
  fk_history       = zeros(nd, M);
  cartesian_errors = zeros(1, M);

  %— main loop ------------------------------------------------------------
  for k = 1:M
    p_k = f_num(qi.');
    err = norm(desired_point - p_k);
    Jk  = J_num(qi.');

    guesses(k, :)       = qi.';
    fk_history(:, k)    = p_k;
    cartesian_errors(k) = err;

    if err < opts.CartesianTolerance, break; end

    method = lower(opts.Method);
    if strcmp(method,'auto')
      if nd==nq
        Qinvg = inv(Jk);
      else
        Qinvg = pinv(Jk);
      end
    else  % dls
      lambda = opts.DLSLambda;
      Qinvg = (Jk.'*Jk + lambda^2*eye(nq)) \ Jk.';
    end

    dq = Qinvg * (desired_point - p_k);
    if norm(dq) < opts.MinJointIncrement, break; end

    if strcmp(method,'auto') && nd==nq && abs(det(Jk))<opts.MaxClosenessSingularity
      warning('Near singularity. Stopping.'); break;
    end

    qi = qi + dq;
  end

  it = k;
  guesses          = guesses(1:it, :);
  fk_history       = fk_history(:, 1:it);
  cartesian_errors = cartesian_errors(1:it);
  q_out = qi;

  if opts.Verbose
    fprintf('\n--- NewtonMethod full history ---\n');
    for i = 1:it
      fprintf('Iter %3d: q = [%s],  f(q) = [%s],  err = %.3e\n', ...
        i, ...
        num2str(guesses(i,:),       '%.6f '), ...
        num2str(fk_history(:,i).', '%.6f '), ...
        cartesian_errors(i) );
    end
    fprintf('----------------------------------\n\n');
  end
end
