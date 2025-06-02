function [q_tau, coefficients, vel_max_time, vel_max_val, acc_max_time, acc_max_val] = quintic_compute_joint_trajectories(q_in, q_fin, T, p_dot_0, p_dot_1, a_in, a_f, link_len, direct_kinematics, print_info)
        % Compute joint trajectories for a multi-joint robot
        %
        % Inputs:
        %   q_in: Initial joint angles (vector)
        %   q_fin: Final joint angles (vector)
        %   T: Total time
        %   p_dot_0: Initial end effector velocity
        %   p_dot_1: Final end effector velocity
        %   a_in: Initial acceleration
        %   a_f: Final acceleration
        %   link_len: Link lengths (vector)
        %   direct_kinematics: Symbolic expression of direct kinematics
        %   print_info : Boolean to print detailed information
        %
        % Outputs:
        %   q_tau: Cell array of symbolic expressions for joint trajectories
        %   coefficients: Matrix of trajectory coefficients (one row per joint)
            % vel_max_time: (istante di tempo ove la velocità raggiunge il
            % massimo
            % vel_max_val: valore di velocità massimo 
            % acc_max_time:; istante di tempo ove l'accelerazione raggiunge il massimo 
            % acc_max_val: valore di accelerazione massimo


        %% esempi di uso: 
        %REST-TO-REST (su posizione) 
        % A = [1;1;1]; 
        % B = [-1;5;0]; 
        % 
        % 
        % [q_tau, coefficients, ... altri output] = quintic_compute_joint_trajectories(A, B, 2.5, [0, 0, 0],[0, 0, 0], [0, 0, 0],[0, 0, 0], [], [], true);

        
        % REST-TO-REST (su orientamento) dove i 3 valori sono il prodotto
        % ottenuto da theta*r estratti da inverseRodriguezSolution.
        % Plotterà quindi variazione di orientamento, velocità angolare e
        % accelerazione angolare. 

        % thetain = [0,0,0]; 
        % thetaf = [-1.73,1.73,0.72]; 
        % 
        % 
        % [q_tau, coefficients, ecc...] = quintic_compute_joint_trajectories(thetain, thetaf, 2.5, [0, 0, 0],[0, 0, 0], [0, 0, 0],[0, 0, 0], [], [], true);
        
        % CASO CON DK. ATTENZIONE LA FORMA DEI P_DOT CAMBIA. 
        % syms q1 q2 q3 L real
        % 
        % f_q = [q2*cos(q1)+0.6*cos(q1+q3);
        %     q2*sin(q1)+0.6*sin(q1+q3);
        %     q1+q3]; 
        % 
        % 
        % [q_tau, coefficients, vel_max_time, vel_max_val, acc_max_time, acc_max_val] = quintic_compute_joint_trajectories(A, B, 2.5, [1;2;3],[0;0;0], [0;0;0],[0;0;0], [], f_q, true);
        %% 
    
        num_joints = length(q_in);
        link_lengths = link_len;
        
        % Jacobian
        if isempty(direct_kinematics)
             % salta blocco Jacobiano
             q_dot_0 = p_dot_0;   % li passi tu da input
             q_dot_T = p_dot_1;
        else
    
            % Substitute initial values for joints and keep other variables symbolic
            % ----- inside the Jacobian branch ---------------------------------
            joint_vars = sym('q', [1 num_joints]);
            all_vars   = symvar(direct_kinematics);
            non_joint_vars = setdiff(all_vars, joint_vars);
            
            J = simplify( jacobian(direct_kinematics, joint_vars) );
            
            % ----- sostituzione --> due casi -------------
            if isempty(link_len)                    % l'utente NON passa L
                J_initial = subs(J, joint_vars, q_in(:).');      % solo i giunti
            else
                % controllo di coerenza:
                if numel(non_joint_vars) ~= numel(link_len)
                    error('link_len deve avere %d elementi (uno per %s).', ...
                          numel(non_joint_vars), char(non_joint_vars));
                end
                J_initial = subs(J, [joint_vars, non_joint_vars], ...
                                     [q_in(:).',  link_len(:).']);
            end

            
            % Check if J_initial is square. If yes, inverse, otherwise we must compute the pseudo-inverse
            if size(J_initial, 1) == size(J_initial, 2)
                J_initial_inv = simplify(inv(J_initial));
            else
                J_initial_inv = simplify(pinv(J_initial));
            end
            disp(J_initial_inv)
    
            % Velocity of configurations
            q_dot_0 = J_initial_inv * p_dot_0;
            q_dot_T = J_initial_inv * p_dot_1;

        end
    
        % Initialize outputs
        q_tau = cell(1, num_joints);
        q_tau_dot = cell(1, num_joints);
        q_tau_dot_dot = cell(1, num_joints);
        coefficients = zeros(num_joints, 6);
    
        % Compute trajectories for each joint
        vel_max_time = zeros(1,num_joints);
        vel_max_val  = zeros(1,num_joints);
        acc_max_time = zeros(1,num_joints);
        acc_max_val  = zeros(1,num_joints);
        
        for i = 1:num_joints
            fprintf('-- Computations for Joint %d:\n', i);
            [q, c] = quintic_poly_double_norm_compute_coeff(q_in(i), q_fin(i), q_dot_0(i), q_dot_T(i), a_in(i),    a_f(i), T, print_info);
            coefficients(i, :) = [c(6), c(5), c(4), c(3), c(2), c(1)];
            q_tau{i} = q;
            q_tau_dot{i} = diff(q, 'tau') / T;
            q_tau_dot_dot{i} = diff(q_tau_dot{i}, 'tau') / T;

            % Derivate in tau
            q_tau_dot_tau  = diff(q, 'tau');          % dq/dτ
            q_tau_dd_tau   = diff(q_tau_dot_tau,'tau');   % d²q/dτ²  (≅ accel. in τ)
            q_tau_ddd_tau  = diff(q_tau_dd_tau,'tau');    % d³q/dτ³  (≅ jerk in τ)
            
            % ------------------- VELOCITY MAX -------------------
            syms tau real
            sol_vel = solve(q_tau_dd_tau == 0, tau);
            sol_vel = double(sol_vel);                       % da sym a double
            sol_vel = sol_vel( imag(sol_vel)==0 & ...        % scarta complessi
                               sol_vel>=0 & sol_vel<=1 );    % dentro [0,1]

            
            cand_tau_vel = [0; 1; sol_vel(:)];           % <-- aggiungi i bordi
            cand_vel = abs( subs(q_tau_dot_tau, tau, cand_tau_vel) ) / T;
            
            [vel_max_val, idx_v] = max( vpa(cand_vel,6) );
            vel_max_tau  = cand_tau_vel(idx_v);
            vel_max_time = vel_max_tau * T;

            
            % ---------------- ACCELERATION MAX ------------------
            sol_acc = solve(q_tau_ddd_tau == 0, tau);
            sol_acc = double(sol_acc);
            sol_acc = sol_acc( imag(sol_acc)==0 & ...
                               sol_acc>=0 & sol_acc<=1 );
            
            cand_tau_acc = [0; 1; sol_acc(:)];           % idem
            cand_acc = abs( subs(q_tau_dd_tau, tau, cand_tau_acc) ) / T^2;
            
            [acc_max_val, idx_a] = max( vpa(cand_acc,6) );
            acc_max_tau  = cand_tau_acc(idx_a);
            acc_max_time = acc_max_tau * T;

            vel_max_time(i) = vel_max_time;
            vel_max_val(i)  = vel_max_val;
            acc_max_time(i) = acc_max_time;
            acc_max_val(i)  = acc_max_val;


            
            % Display results
            if(print_info == false)
                disp(q_tau{i});
                fprintf('\n');
            end
        end
    
        % Plot trajectories
        disp('ATTENZIONE: GLI ISTANTI DI TEMPO MASSIMI PER VELOCITA E ACCELERAZIONE SONO RESTITUITI COME  t=tau*T')
        plot_trajectories(q_tau, q_tau_dot, q_tau_dot_dot, T);
    end
        
    function plot_trajectories(q_tau, q_tau_dot, q_tau_dot_dot, T)
        num_joints = length(q_tau);
        syms t real
        tau = t/T;
        
        % Position
        figure;
        hold on;
        for i = 1:num_joints
            fplot(subs(q_tau{i}, 'tau', tau), [0, T]);
        end
        title('Joint Positions');
        xlabel('Time (s)');
        ylabel('Position (rad)');
        legend(arrayfun(@(x) sprintf('q%d', x), 1:num_joints, 'UniformOutput', false));
        hold off;
    
        % Velocity
        figure;
        hold on;
        for i = 1:num_joints
            fplot(subs(q_tau_dot{i}, 'tau', tau), [0, T]);
        end
        title('Joint Velocities');
        xlabel('Time (s)');
        ylabel('Velocity (rad/s)');
        h = legend(arrayfun(@(x) sprintf('q%d_dot', x), 1:num_joints, 'UniformOutput', false));
        set(h, 'Interpreter', 'none');
        hold off;
    
        % Acceleration
        figure;
        hold on;
        for i = 1:num_joints
            fplot(subs(q_tau_dot_dot{i}, 'tau', tau), [0, T]);
        end
        title('Joint Accelerations');
        xlabel('Time (s)');
        ylabel('Acceleration (rad/s^2)');
        h = legend(arrayfun(@(x) sprintf('q%d_dot_dot', x), 1:num_joints, 'UniformOutput', false));
        set(h, 'Interpreter', 'none');
        hold off;
    end