% This code will perform an optimization on v_min, v_max, and lift_threshold 
% to maximize the final horizontal distance (final x position) of the projectile.

% --- SETUP: Parallel pool ---
% Uncomment the following line if you don't have a parallel pool open
% parpool; 

% Load regression coefficient tables
cd_regression_coefficients = 'C_d_regression_coefficients.csv';
cl_regression_coefficients = 'C_l_regression_coefficients.csv';
cd_coeffTable = readtable(cd_regression_coefficients);
cl_coeffTable = readtable(cl_regression_coefficients);

% Create symbolic variables for Mach and AoA
syms Mach AoA

% Construct the symbolic polynomials
cd_poly = construct_symbolic_poly(cd_coeffTable, Mach, AoA);
cl_poly = construct_symbolic_poly(cl_coeffTable, Mach, AoA);

% Define the objective function for optimization
% params = [v_min, v_max, lift_threshold]
objFun = @(params) simulate_projectile(params(1), params(2), params(3), cd_poly, cl_poly);

% Set bounds for v_min, v_max, and lift_threshold
% You may adjust these bounds based on what makes sense for your scenario.
lb = [-200, 0, 5];    % Lower bounds
ub = [-50, 100, 200];  % Upper bounds (example: lift_threshold up to 200 N above m*g)

% Adjust optimization options to change stopping criteria
options = optimoptions('fmincon', 'Display', 'iter', ...
    'UseParallel', false, ...    % We'll handle parallelization manually
    'OptimalityTolerance', 1e-6, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 10, ...
    'MaxFunctionEvaluations', 5000);

% Number of random initial guesses
numStarts = 3;

% Generate random initial guesses within bounds
initial_guesses = lb + rand(numStarts, length(lb)) .* (ub - lb);

% Preallocate arrays to hold results
optimal_params_array = zeros(numStarts, length(lb));
fval_array = zeros(numStarts, 1);

% Run optimizations in parallel
parfor i = 1:numStarts
    % Random initial guess
    initial_guess = initial_guesses(i, :);
    
    % Create optimization problem
    problem = createOptimProblem('fmincon', 'objective', objFun, ...
        'x0', initial_guess, 'lb', lb, 'ub', ub, 'options', options);
    
    % Run optimization
    [x, fval] = fmincon(problem);
    
    % Store results
    optimal_params_array(i, :) = x;
    fval_array(i) = fval;
end

% Find the best result
[~, bestIdx] = min(fval_array);
optimal_params = optimal_params_array(bestIdx, :);
fval = fval_array(bestIdx);

v_min_opt = optimal_params(1);
v_max_opt = optimal_params(2);
lift_threshold_opt = optimal_params(3);

fprintf('Optimized v_min: %.2f\n', v_min_opt);
fprintf('Optimized v_max: %.2f\n', v_max_opt);
fprintf('Optimized lift_threshold: %.2f\n', lift_threshold_opt);
fprintf('Maximum final X position achieved: %.2f km\n', -fval/1000);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = simulate_projectile(v_min, v_max, lift_threshold, cd_poly, cl_poly)
    % This function runs the projectile simulation with given v_min, v_max, and 
    % lift_threshold, and returns the negative final x position.
    %
    % Inputs:
    %  - v_min: velocity threshold for trigger condition
    %  - v_max: velocity threshold for resetting trigger
    %  - lift_threshold: additional lift force threshold above m*g
    %  - cd_poly, cl_poly: symbolic polynomials for drag and lift
    %
    % Output:
    %  - cost: negative of final x position to be minimized by optimizer

    aoa_min = 0; 
    aoa_max = 15;

    M = 10; % Initial Mach number
    m = 25; % mass (kg)
    y = 40000;
    x = 0;
    [T, a, P, rho] = atmosisa(y, 'extended', true);

    A = 2; % Surface area
    back_area = 0.028; %m^2

    vx = M*a;
    vy = 0;

    g = 9.81;
    dt = 1;

    trigger = false;

    while y > 0
        [T, a, P, rho] = atmosisa(y, 'extended', true);

        v = sqrt(vx^2 + vy^2);
        M = v/a;

        L_fun = @(aoa) double(subs(cl_poly, {sym('Mach'), sym('AoA')}, [M, aoa]));
        D_fun = @(aoa) double(subs(cd_poly, {sym('Mach'), sym('AoA')}, [M, aoa]));

        if vy < v_min && vy > v_min*1.2
            trigger = true;
        end

        if vy > v_max
            trigger = false;
        end

        % L_min: lift required to exceed gravity by some threshold
        L_min = m*g + lift_threshold;

        if trigger
            % We want to maximize L/D but ensure L >= L_min.
            % obj_fun: maximize L/D => minimize -L/D.
            obj_fun = @(aoa) -L_fun(aoa)/D_fun(aoa);

            % Nonlinear constraint ensures L(aoa) >= L_min
            nonlcon = @(aoa) deal([], L_fun(aoa) - L_min);

            aoa0 = (aoa_min + aoa_max)/2;
            opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
            [aoa, ~] = fmincon(obj_fun, aoa0, [], [], [], [], aoa_min, aoa_max, nonlcon, opts);
        else
            % If not triggered, just maximize L/D (no constraint)
            obj_fun = @(aoa) -(L_fun(aoa)/D_fun(aoa));
            [aoa, ~] = fminbnd(obj_fun, aoa_min, aoa_max);
        end

        L = L_fun(aoa);
        D = D_fun(aoa);

        Lx = -L * vy / v;
        Ly = L * vx / v;

        Dx = -D * vx/v;
        Dy = -D * vy/v;

        Fx = Dx+Lx;
        Fy = Ly + Dy - m*g;

        vx = vx + (Fx/m)*dt;
        vy = vy + (Fy/m)*dt;

        x = x + vx*dt;
        y = y + vy*dt;
    end

    % cost is negative of final x
    cost = -x;
end

function poly = construct_symbolic_poly(coeffTable, Mach, AoA)
    % Initialize the polynomial
    poly = 0;
    
    % Loop through each row in the table
    for i = 1:height(coeffTable)
        feature = coeffTable.Feature{i};
        coefficient = coeffTable.Coefficient(i);
        
        % Start with the coefficient
        term = coefficient;
        
        % Check if the feature is not the intercept
        if ~strcmp(feature, 'Intercept')
            % Split the feature into parts (e.g., 'Mach', 'Mach^2', 'Mach AoA')
            terms = split(feature, ' ');
            for j = 1:length(terms)
                if contains(terms{j}, '^')
                    % Extract the variable and its exponent (e.g., 'Mach^2')
                    [var, exponent] = strtok(terms{j}, '^');
                    exponent = str2double(extractAfter(exponent, '^'));
                else
                    % Default exponent is 1 for terms like 'Mach' or 'AoA'
                    var = terms{j};
                    exponent = 1;
                end
                
                % Multiply the term by the corresponding variable raised to its exponent
                if strcmp(var, 'Mach')
                    term = term * Mach^exponent;
                elseif strcmp(var, 'AoA')
                    term = term * AoA^exponent;
                end
            end
        end
        
        % Add the term to the polynomial
        poly = poly + term;
    end
end
