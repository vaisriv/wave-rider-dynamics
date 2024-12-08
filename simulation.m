% Clear history
clear;
clc;

% Set LaTeX as Interpreter for text labels
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

% Trial Directories
trial_dirs = [
    './trial1/';
    './trial2/';
    './trial3/';
    './trial4/';
];
trial_num = 1;

% Read the CSV files into tables
cd_coeffTable = readtable([trial_dirs(trial_num, :) 'C_d_Regression.csv']);
cl_coeffTable = readtable([trial_dirs(trial_num, :) 'C_l_Regression.csv']);

% Create symbolic variables for Mach and AoA
syms Mach AoA

% Generate symbolic polynomial for C_d, C_l, and C_l/C_d
cd_poly = construct_symbolic_poly(cd_coeffTable, Mach, AoA);
cl_poly = construct_symbolic_poly(cl_coeffTable, Mach, AoA);
cl_cd_poly = cl_poly/cd_poly;

% Generate C_l/C_d for all AoA
aoa_linsp = linspace(0, 15, 31);
cl_cd_mach_poly = subs(cl_cd_poly, AoA, aoa_linsp);

% Intial Conditions
M = 10; % Initial mach number
m = 25; %25 kg
y = 40000;
x = 0;
aoa = 0;
[T, a, P, rho] = atmosisa(y, 'extended', true);

A = 2; % Surface area
back_areas = [
    0.028;
    0.030;
    0.030;
    0.0298;
]; %m^2

vx = M*a;
vy = 0;

g = 9.81;
dt = 1;

trajectory_plot = figure;
hold on;
xlabel('X Position');

yyaxis left; % put trajectory on left
traj_plt = plot(NaN, NaN, 'b-'); % Empty plot for trajectory
ylabel('Y Position');
ylim([0 4.5e4]);

yyaxis right; % put mach on right
mach_plt = plot(NaN, NaN, 'r-'); % Empty plot for mach number
ylabel('Mach Number');
ylim([0 11]);

title('Waverider Trajectory Simulation');
plotHandle = plot(x, y, 'bo');

% initial simulation state
trajectoryX = x;
trajectoryY = y;
M_history = M;
aoa_history = aoa;

last_saved_x = x;

% cutoff height for when mach number is low enough that shock is weak and it is time to glide
glide_cutoff = 2.5e4;

while y > 0

    if y > glide_cutoff
        [T, a, P, rho] = atmosisa(y, 'extended', true);
    else
        [T, a, P, rho] = atmosisa(glide_cutoff, 'extended', true);
    end

    v = sqrt(vx^2 + vy^2);
    M = v/a;

    if x == 0 || y < glide_cutoff
        aoa = 0;
        L = double(subs(cl_poly, [Mach AoA], [M aoa])*rho);
        D = double(subs(cd_poly, [Mach AoA], [M aoa])*rho)+P/M*back_areas(trial_num);
        cl_cd = L/D;
    else
        [aoa, L, D, cl_cd] = getMaxCLCD(aoa_linsp, cl_cd_mach_poly, cl_poly, cd_poly, Mach, AoA, M, rho, P, back_areas(trial_num));
    end

    Lx = -L * vy / v;
    Ly = L * vx / v;

    Dx = - D * vx/v;
    Dy = - D * vy/v;

    Fx = Dx+Lx;
    Fy = Ly + Dy -m*g;

    vx = vx + (Fx/m)*dt;
    vy = vy + (Fy/m)*dt;

    x = x + vx*dt;
    y = y + vy*dt;

    % Store AoA and Mach in history arrays
    aoa_history = [aoa_history, aoa];
    M_history = [M_history, M];

    trajectoryX = [trajectoryX, x];
    trajectoryY = [trajectoryY, y];
    last_saved_x = x;
    
    % Update the plot data
    set(traj_plt, 'XData', trajectoryX, 'YData', trajectoryY);
    set(mach_plt, 'XData', trajectoryX, 'YData', M_history);
    drawnow;

    if abs(x - last_saved_x) > 1000
        % Print information
        fprintf('AoA: %.2f degrees\n', aoa);
        fprintf('Mach Number: %.2f\n', M);
        fprintf('Lift: %.4f\n', L);
        fprintf('Drag: %.4f\n', D);
        fprintf('$\frac{C_L}{C_D}$: %.4f\n', cl_cd);
        fprintf('$L_x$: %.2f N, $L_y$: %.2f N\n', Lx, Ly);
        fprintf('Air Density ($\rho$): %.5f kg/m^3\n', rho);
        fprintf('Velocity (v): %.2f m/s\n', v);
        fprintf('$F_x$: %.2f N, $F_y$: %.2f N\n', Fx, Fy);
        fprintf('$v_x$: %.2f m/s, $v_y$: %.2f m/s\n', vx, vy);
        fprintf('----------------------------------------\n');
    end
end

fprintf('Final X Position: %.2f kilometers\n', x/1000);
set(mach_plt, 'XData', trajectoryX, 'YData', M_history);
drawnow;

saveas(trajectory_plot, [trial_dirs(trial_num, :) 'trajectory.png']);
hold off;

%{
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
%}

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

function [aoa, L, D, cl_cd] = getMaxCLCD(aoa_linsp, cl_cd_mach_poly, cl_poly, cd_poly, Mach, AoA, M, rho, P, back_area)
    [cl_cd, i] = max(subs(cl_cd_mach_poly, Mach, M));
    cl_cd = double(cl_cd);
    aoa = aoa_linsp(i);
    L = double(subs(cl_poly, [Mach AoA], [M aoa])*rho);
    D = double(subs(cd_poly, [Mach AoA], [M aoa])*rho)+P/M*back_area;
end