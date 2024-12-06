% Specify the path to your CSV files
cd_regression_coefficients = 'C_d_regression_coefficients.csv';
cl_regression_coefficients = 'C_l_regression_coefficients.csv';

% Read the CSV files into tables
cd_coeffTable = readtable(cd_regression_coefficients);
cl_coeffTable = readtable(cl_regression_coefficients);

% Create symbolic variables for Mach and AoA
syms Mach AoA

% Generate symbolic polynomial for C_d
cd_poly = construct_symbolic_poly(cd_coeffTable, Mach, AoA);
fprintf('Symbolic Polynomial for C_d:\n');
disp(cd_poly);

% Generate symbolic polynomial for C_l
cl_poly = construct_symbolic_poly(cl_coeffTable, Mach, AoA);
fprintf('Symbolic Polynomial for C_l:\n');
disp(cl_poly);

% Define new input data (Mach, AoA)
new_data = [0.8, 5.0]; % Example: Mach 0.8, AoA 5 degrees

% Calculate values of the polynomials at given inputs
predicted_cd = double(subs(cd_poly, {Mach, AoA}, new_data));
predicted_cl = double(subs(cl_poly, {Mach, AoA}, new_data));

% This is the simulation file for the 2D projectile motion simulation

aoa_min = 0; 
aoa_max = 15;

M = 10; % Initial mach number
m = 25; %25 kg
y = 40000;
x = 0;
[T, a, P, rho] = atmosisa(y, 'extended', true);

A = 2; % Surface area
back_area = 0.028; %m^2

vx = M*a;
vy = 0;

g = 9.81;
dt = 1;

figure;
hold on;
xlabel('X Position');
ylabel('Y Position');
title('Waverider Trajectory Simulation');
plotHandle = plot(x, y, 'bo');
trajectoryX = x;
trajectoryY = y;

% Initialize the plot handle outside the loop
h = plot(NaN, NaN, 'b-'); % Create an empty plot

last_saved_x = x;
min_v = -100;
maax_v = 40;

trigger = false;

aoa = 0;

lift_threshold = 20;

% Initialize arrays to store AoA and Mach during simulation
aoa_history = [];
M_history = [];

while y > 0
    [T, a, P, rho] = atmosisa(y, 'extended', true);

    v = sqrt(vx^2 + vy^2);
    M = v/a;

    L = @(aoa) double(subs(cl_poly, {Mach, AoA}, [M, aoa]));
    D = @(aoa) double(subs(cd_poly, {Mach, AoA}, [M, aoa]));

    if vy < min_v && vy > min_v*1.2
        trigger = true;
    end
    
    if vy > maax_v
        trigger = false;
    end

    L_min = m*g+lift_threshold;

    if trigger
        % Define the objective function (maximize L/D, so minimize -L/D)
        obj_fun = @(aoa) -L(aoa) / D(aoa);

        % Define the nonlinear constraint
        % Ensure L(aoa) >= L_min
        nonlcon = @(aoa) deal([], L(aoa) - L_min); % No equality constraint, inequality only

        % Initial guess for aoa
        aoa0 = (aoa_min + aoa_max) / 2; 

        % Set options for fmincon
        options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp'); 

        % Call fmincon to optimize
        [aoa, fval] = fmincon(obj_fun, aoa0, [], [], [], [], aoa_min, aoa_max, nonlcon, options);
    else
        obj_fun = @(aoa) -(L(aoa) / D(aoa));
        [aoa, fval] = fminbnd(obj_fun, aoa_min, aoa_max);
    end

    % Recompute L and D with the chosen AoA
    L_val = double(subs(cl_poly, {Mach, AoA}, [M, aoa]));
    D_val = double(subs(cd_poly, {Mach, AoA}, [M, aoa]));

    Lx = -L_val * vy / v;
    Ly = L_val * vx / v;

    Dx = - D_val * vx/v;
    Dy = - D_val * vy/v;

    Fx = Dx+Lx;
    Fy = Ly + Dy -m*g;

    vx = vx + (Fx/m)*dt;
    vy = vy + (Fy/m)*dt;

    x = x + vx*dt;
    y = y + vy*dt;

    % Store AoA and Mach in history arrays
    aoa_history = [aoa_history, aoa];
    M_history = [M_history, M];

    if abs(x - last_saved_x) > 1000
        trajectoryX = [trajectoryX, x];
        trajectoryY = [trajectoryY, y];
        last_saved_x = x;

        % Print information
        fprintf('Optimal AoA: %.2f degrees\n', aoa);
        fprintf('Mach Number: %.2f\n', M);
        fprintf('Lift: %.4f\n', L_val);
        fprintf('Drag: %.4f\n', D_val);
        fprintf('Lift (L): %.2f N\n', L_val);
        fprintf('Drag (D): %.2f N\n', D_val);
        fprintf('Lx: %.2f N, Ly: %.2f N\n', Lx, Ly);
        fprintf('Air Density (rho): %.5f kg/m^3\n', rho);
        fprintf('Velocity (v): %.2f m/s\n', v);
        fprintf('Fx: %.2f N, Fy: %.2f N\n', Fx, Fy);
        fprintf('vx: %.2f m/s, vy: %.2f m/s\n', vx, vy);
        fprintf('Trigger: %d\n', trigger);
        fprintf('----------------------------------------\n');

        % Update the plot data
        set(h, 'XData', trajectoryX, 'YData', trajectoryY);
        drawnow;
    end
end

fprintf('Final X Position: %.2f kilometers\n', x/1000);

hold off;

% After the simulation, plot the distributions of AoA and Mach
figure;
histogram(aoa_history);
xlabel('Angle of Attack (degrees)');
ylabel('Frequency');
title('Distribution of AoA During Flight');

figure;
histogram(M_history);
xlabel('Mach Number');
ylabel('Frequency');
title('Distribution of Mach Number During Flight');

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
