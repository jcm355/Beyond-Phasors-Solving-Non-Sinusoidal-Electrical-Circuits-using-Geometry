%==========================================================================
% COMPUTATIONAL BENCHMARKING PROGRAM: 4-CASE (Tesis / Subplots)
% Version: 24.4-tesis-subplots
%==========================================================================
clear all; clc;

% --- SECTION 1: GA-FuL INITIALIZATION (Only for PARSING) ---
try
    addpath('v2025-07-16_Libreria_Egipcio');
    library_folder = fullfile('v2025-07-16_Libreria_Egipcio', 'Release');
    gapotAssemblyPath = fullfile(pwd, library_folder);
    addpath(gapotAssemblyPath);
    
    global gafulAssembly ga
    dotnetenv("framework");
    if (isempty(gafulAssembly))
        gafulAssembly = NET.addAssembly(fullfile(gapotAssemblyPath, 'GeometricAlgebraFulcrumLib.Matlab.exe'));
    end
    fprintf('GA-FuL Library loaded (only for input parsing).\n');
catch ME
    error('Could not load GA-FuL library. Check the path: %s', ME.message);
end

disp('======================================================');
disp('   COMPUTATIONAL BENCHMARKING: 4-CASE (v24.4-tesis-subplots)   ');
disp(' (GA-Real (Native) vs. Phasors (Native)) ');
disp('======================================================');

% --- User Inputs for Benchmarking ---
N_max = input('Enter the maximum number of harmonics to test (N_max): ');
N_step = input('Enter the harmonic step (T) (e.g., 10 to test 1, 11, 21...): ');
if isempty(N_step) || N_step < 1
    N_step = 1;
end
num_mediciones = input('Enter the number of measurements per point (b) (e.g., 5): ');
if isempty(num_mediciones) || num_mediciones < 1
    num_mediciones = 1;
end

% --- Common Circuit Parameters ---
disp('--- Enter Common Circuit Parameters ---');
R_val = input('Resistance R (Ohms): ');
L_val = input('Inductance L (Henries): ');
C_val = input('Capacitance C (Farads): ');
omega = input('Fundamental angular frequency w (rad/s): ');

% --- Common GA-FuL Setup ---
dimension_max = N_max * 2;
ga = gafulGetProcessor(dimension_max, 0); 
disp('ATTENTION: Use GA-FuL format for input (0-based indices)!');

% --- Input Vector 1: VOLTAGE (for Cases 1a, 2b) ---
voltage_input_text = input('Enter the 1st H. VOLTAGE vector (e.g., 1<0>): ', 's');
% --- Input Vector 2: CURRENT (for Cases 1b, 2a) ---
current_input_text = input('Enter the 1st H. CURRENT vector (e.g., 1.5<0>): ', 's');

decay_factor = input('Enter the decay factor (e.g., 0.1 for 10%): ');
iteraciones_bulk_nativos = 10000; 
fprintf('A batch of %d iterations will be used for NATIVE (Phasor, GA-Real).\n', iteraciones_bulk_nativos);
fprintf('%d timeit measurements will be taken for each point N and averaged.\n', num_mediciones);

% --- Harmonics Axis (N_axis) ---
N_axis = 1:N_step:N_max;
if N_axis(end) ~= N_max 
    N_axis = [N_axis, N_max];
end
num_steps = length(N_axis);

% --- Pre-allocate storage for results ---
% Columns: 1=(1a), 2=(1b), 3=(2a), 4=(2b)
tiempos_ga_pure = zeros(num_steps, 4);
tiempos_clasico_pure = zeros(num_steps, 4);
std_ga_pure = zeros(num_steps, 4);
std_clasico_pure = zeros(num_steps, 4);

% --- Generation of FULL VOLTAGE vector (for 1a, 2b) ---
vector_fundamental_u = ga.Parse(voltage_input_text);
[coef_fundamental_u] = parsear_vector_ga_ful(vector_fundamental_u, 2);
A1_u = coef_fundamental_u(1); 
B1_u = coef_fundamental_u(2); 
coeficientes_totales_u = zeros(1, dimension_max);
scale_factor_u = 1.0;
for h = 1:N_max
    coeficientes_totales_u(2*h - 1) = A1_u * scale_factor_u;
    coeficientes_totales_u(2*h) = B1_u * scale_factor_u;
    scale_factor_u = scale_factor_u * (1 - decay_factor);
end
coef_cos_u_full = coeficientes_totales_u(1:2:end);
coef_sen_u_full = coeficientes_totales_u(2:2:end);
fprintf('Voltage input vector generated for N=%d.\n', N_max);

% --- Generation of FULL CURRENT vector (for 1b, 2a) ---
vector_fundamental_i = ga.Parse(current_input_text);
[coef_fundamental_i] = parsear_vector_ga_ful(vector_fundamental_i, 2);
A1_i = coef_fundamental_i(1); 
B1_i = coef_fundamental_i(2); 
coeficientes_totales_i = zeros(1, dimension_max);
scale_factor_i = 1.0;
for h = 1:N_max
    coeficientes_totales_i(2*h - 1) = A1_i * scale_factor_i;
    coeficientes_totales_i(2*h) = B1_i * scale_factor_i;
    scale_factor_i = scale_factor_i * (1 - decay_factor);
end
coef_cos_i_full = coeficientes_totales_i(1:2:end);
coef_sen_i_full = coeficientes_totales_i(2:2:end);
fprintf('Current input vector generated for N=%d.\n', N_max);

% --- Pre-calculate Reactance/Susceptance vectors ---
k_indices_max = 1:N_max;
Xk_vector_max = k_indices_max .* omega .* L_val - 1./(k_indices_max .* omega .* C_val);
Xk_vector_max(isnan(Xk_vector_max)) = 0;

G_val = 1 / R_val;
Bk_vector_max = k_indices_max .* omega .* C_val - 1./(k_indices_max .* omega .* L_val);
Bk_vector_max(isnan(Bk_vector_max)) = 0;

disp('======================================================');
disp('STARTING BENCHMARK SWEEP...');
disp('======================================================');

%==========================================================================
% --- CASE 1a: SERIES RLC (u -> i) ---
%==========================================================================
disp('--- Running Case 1a: Series RLC (u -> i) ---');
fprintf('Starting sweep in %d steps (from N=1 to N=%d)...\n', num_steps, N_max);
for i = 1:num_steps
    n = N_axis(i); 
    
    % Get coefficients for this N
    coef_cos_u = coef_cos_u_full(1:n);
    coef_sen_u = coef_sen_u_full(1:n);
    Xk_vector = Xk_vector_max(1:n);
    
    temp_ga = zeros(1, num_mediciones);
    temp_clasico = zeros(1, num_mediciones);
    
    fprintf('  Step %d/%d (N = %d)... measuring %d times:\n', i, num_steps, n, num_mediciones);
    for j = 1:num_mediciones
        % 1. Time GA-Real
        f_ga = @() loop_ga_1a(n, coef_cos_u, coef_sen_u, Xk_vector, R_val, iteraciones_bulk_nativos);
        temp_ga(j) = timeit(f_ga) / iteraciones_bulk_nativos;
        
        % 2. Time Phasor
        f_clasico = @() loop_fasor_1a(n, coef_cos_u, coef_sen_u, Xk_vector, R_val, iteraciones_bulk_nativos);
        temp_clasico(j) = timeit(f_clasico) / iteraciones_bulk_nativos;
    end
    % --- Calculate and store average and std dev ---
    tiempos_ga_pure(i, 1) = mean(temp_ga);
    tiempos_clasico_pure(i, 1) = mean(temp_clasico);
    std_ga_pure(i, 1) = std(temp_ga);
    std_clasico_pure(i, 1) = std(temp_clasico);
    
    fprintf('    GA-Real: %.3e (std %.2e) | Phasor: %.3e (std %.2e)\n', ...
        tiempos_ga_pure(i, 1), std_ga_pure(i, 1), ...
        tiempos_clasico_pure(i, 1), std_clasico_pure(i, 1));
end 

%==========================================================================
% --- CASE 1b: SERIES RLC (i -> u) ---
%==========================================================================
disp('--- Running Case 1b: Series RLC (i -> u) ---');
fprintf('Starting sweep in %d steps (from N=1 to N=%d)...\n', num_steps, N_max);
for i = 1:num_steps
    n = N_axis(i);
    
    % Get coefficients for this N
    coef_cos_i = coef_cos_i_full(1:n);
    coef_sen_i = coef_sen_i_full(1:n);
    Xk_vector = Xk_vector_max(1:n);
    
    temp_ga = zeros(1, num_mediciones);
    temp_clasico = zeros(1, num_mediciones);
    
    fprintf('  Step %d/%d (N = %d)... measuring %d times:\n', i, num_steps, n, num_mediciones);
    for j = 1:num_mediciones
        f_ga = @() loop_ga_1b(n, coef_cos_i, coef_sen_i, Xk_vector, R_val, iteraciones_bulk_nativos);
        temp_ga(j) = timeit(f_ga) / iteraciones_bulk_nativos;
        
        f_clasico = @() loop_fasor_1b(n, coef_cos_i, coef_sen_i, Xk_vector, R_val, iteraciones_bulk_nativos);
        temp_clasico(j) = timeit(f_clasico) / iteraciones_bulk_nativos;
    end
    
    % --- Calculate and store average and std dev ---
    tiempos_ga_pure(i, 2) = mean(temp_ga);
    tiempos_clasico_pure(i, 2) = mean(temp_clasico);
    std_ga_pure(i, 2) = std(temp_ga);
    std_clasico_pure(i, 2) = std(temp_clasico);
    
    fprintf('    GA-Real: %.3e (std %.2e) | Phasor: %.3e (std %.2e)\n', ...
        tiempos_ga_pure(i, 2), std_ga_pure(i, 2), ...
        tiempos_clasico_pure(i, 2), std_clasico_pure(i, 2));
end 

%==========================================================================
% --- CASE 2a: PARALLEL RLC (i -> u) ---
%==========================================================================
disp('--- Running Case 2a: Parallel RLC (i -> u) ---');
fprintf('Starting sweep in %d steps (from N=1 to N=%d)...\n', num_steps, N_max);
for i = 1:num_steps
    n = N_axis(i);
    
    % Get coefficients for this N
    coef_cos_i = coef_cos_i_full(1:n);
    coef_sen_i = coef_sen_i_full(1:n);
    Bk_vector = Bk_vector_max(1:n);
    
    temp_ga = zeros(1, num_mediciones);
    temp_clasico = zeros(1, num_mediciones);
    
    fprintf('  Step %d/%d (N = %d)... measuring %d times:\n', i, num_steps, n, num_mediciones);
    for j = 1:num_mediciones
        f_ga = @() loop_ga_2a(n, coef_cos_i, coef_sen_i, Bk_vector, G_val, iteraciones_bulk_nativos);
        temp_ga(j) = timeit(f_ga) / iteraciones_bulk_nativos;
        
        f_clasico = @() loop_fasor_2a(n, coef_cos_i, coef_sen_i, Bk_vector, G_val, iteraciones_bulk_nativos);
        temp_clasico(j) = timeit(f_clasico) / iteraciones_bulk_nativos;
    end
    
    % --- Calculate and store average and std dev ---
    tiempos_ga_pure(i, 3) = mean(temp_ga);
    tiempos_clasico_pure(i, 3) = mean(temp_clasico);
    std_ga_pure(i, 3) = std(temp_ga);
    std_clasico_pure(i, 3) = std(temp_clasico);
    
    fprintf('    GA-Real: %.3e (std %.2e) | Phasor: %.3e (std %.2e)\n', ...
        tiempos_ga_pure(i, 3), std_ga_pure(i, 3), ...
        tiempos_clasico_pure(i, 3), std_clasico_pure(i, 3));
end 

%==========================================================================
% --- CASE 2b: PARALLEL RLC (u -> i) ---
%==========================================================================
disp('--- Running Case 2b: Parallel RLC (u -> i) ---');
fprintf('Starting sweep in %d steps (from N=1 to N=%d)...\n', num_steps, N_max);
for i = 1:num_steps
    n = N_axis(i);
    
    % Get coefficients for this N
    coef_cos_u = coef_cos_u_full(1:n);
    coef_sen_u = coef_sen_u_full(1:n);
    Bk_vector = Bk_vector_max(1:n);
    
    temp_ga = zeros(1, num_mediciones);
    temp_clasico = zeros(1, num_mediciones);
    
    fprintf('  Step %d/%d (N = %d)... measuring %d times:\n', i, num_steps, n, num_mediciones);
    for j = 1:num_mediciones
        f_ga = @() loop_ga_2b(n, coef_cos_u, coef_sen_u, Bk_vector, G_val, iteraciones_bulk_nativos);
        temp_ga(j) = timeit(f_ga) / iteraciones_bulk_nativos;
        
        f_clasico = @() loop_fasor_2b(n, coef_cos_u, coef_sen_u, Bk_vector, G_val, iteraciones_bulk_nativos);
        temp_clasico(j) = timeit(f_clasico) / iteraciones_bulk_nativos;
    end
    
    % --- Calculate and store average and std dev ---
    tiempos_ga_pure(i, 4) = mean(temp_ga);
    tiempos_clasico_pure(i, 4) = mean(temp_clasico);
    std_ga_pure(i, 4) = std(temp_ga);
    std_clasico_pure(i, 4) = std(temp_clasico);
    
    fprintf('    GA-Real: %.3e (std %.2e) | Phasor: %.3e (std %.2e)\n', ...
        tiempos_ga_pure(i, 4), std_ga_pure(i, 4), ...
        tiempos_clasico_pure(i, 4), std_clasico_pure(i, 4));
end 

disp('======================================================');
disp('SWEEP COMPLETE. Generating Plots and Tables...');
disp('======================================================');

%==========================================================================
% --- POST-PROCESSING AND PLOTTING SECTION (4-in-1 Subplots) ---
%==========================================================================

% --- Smoothing factor ---
k_mov = 1; 
if N_step == 1
    k_mov = 5; 
end

% --- Create ONE figure for all 4 subplots ---
figure('Name', 'Combined Benchmark Results (4-Case)', 'NumberTitle', 'off', 'WindowState', 'maximized');

% --- Subplot 1: Case 1a ---
ax1 = subplot(2, 2, 1);
t_ga_smooth = fillmissing(movmean(tiempos_ga_pure(:, 1), k_mov), 'linear', 'EndValues','nearest');
t_cl_smooth = fillmissing(movmean(tiempos_clasico_pure(:, 1), k_mov), 'linear', 'EndValues','nearest');

plot(ax1, N_axis, t_ga_smooth, 'r-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'GA-Real');
hold(ax1, 'on');
plot(ax1, N_axis, t_cl_smooth, 'b-x', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Phasors');
hold(ax1, 'off');
title(ax1, 'Case 1a: Series RLC (u -> i)');
ylabel(ax1, 'Time (seconds)');
legend(ax1, 'show', 'Location', 'NorthWest');
grid(ax1, 'on');
set(ax1, 'YScale', 'log'); 

% --- Subplot 2: Case 1b ---
ax2 = subplot(2, 2, 2);
t_ga_smooth = fillmissing(movmean(tiempos_ga_pure(:, 2), k_mov), 'linear', 'EndValues','nearest');
t_cl_smooth = fillmissing(movmean(tiempos_clasico_pure(:, 2), k_mov), 'linear', 'EndValues','nearest');

plot(ax2, N_axis, t_ga_smooth, 'r-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'GA-Real');
hold(ax2, 'on');
plot(ax2, N_axis, t_cl_smooth, 'b-x', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Phasors');
hold(ax2, 'off');
title(ax2, 'Case 1b: Series RLC (i -> u)');
ylabel(ax2, 'Time (seconds)');
legend(ax2, 'show', 'Location', 'NorthWest');
grid(ax2, 'on');
set(ax2, 'YScale', 'log'); 

% --- Subplot 3: Case 2a ---
ax3 = subplot(2, 2, 3);
t_ga_smooth = fillmissing(movmean(tiempos_ga_pure(:, 3), k_mov), 'linear', 'EndValues','nearest');
t_cl_smooth = fillmissing(movmean(tiempos_clasico_pure(:, 3), k_mov), 'linear', 'EndValues','nearest');

plot(ax3, N_axis, t_ga_smooth, 'r-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'GA-Real');
hold(ax3, 'on');
plot(ax3, N_axis, t_cl_smooth, 'b-x', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Phasors');
hold(ax3, 'off');
title(ax3, 'Case 2a: Parallel RLC (i -> u)');
xlabel(ax3, sprintf('Number of Harmonics (N) [Steps T=%d]', N_step));
ylabel(ax3, 'Time (seconds)');
legend(ax3, 'show', 'Location', 'NorthWest');
grid(ax3, 'on');
set(ax3, 'YScale', 'log');

% --- Subplot 4: Case 2b ---
ax4 = subplot(2, 2, 4);
t_ga_smooth = fillmissing(movmean(tiempos_ga_pure(:, 4), k_mov), 'linear', 'EndValues','nearest');
t_cl_smooth = fillmissing(movmean(tiempos_clasico_pure(:, 4), k_mov), 'linear', 'EndValues','nearest');

plot(ax4, N_axis, t_ga_smooth, 'r-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'GA-Real');
hold(ax4, 'on');
plot(ax4, N_axis, t_cl_smooth, 'b-x', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Phasors');
hold(ax4, 'off');
title(ax4, 'Case 2b: Parallel RLC (u -> i)');
xlabel(ax4, sprintf('Number of Harmonics (N) [Steps T=%d]', N_step));
ylabel(ax4, 'Time (seconds)');
legend(ax4, 'show', 'Location', 'NorthWest');
grid(ax4, 'on');
set(ax4, 'YScale', 'log');

% --- Add a main title to the entire figure ---
sgtitle('Pure Calculation Comparison (Algorithm Efficiency)', 'FontSize', 14, 'FontWeight', 'bold');

% --- Save the combined figure ---
saveas(gcf, 'plot_combined_4_cases.png');
fprintf('Combined 4-in-1 Plot saved as: plot_combined_4_cases.png\n');


%==========================================================================
% --- RESULTS TABLE GENERATION ---
%==========================================================================

% --- Table 1: Case 1a ---
T_1a = create_results_table(N_axis', ...
    tiempos_clasico_pure(:, 1), std_clasico_pure(:, 1), ...
    tiempos_ga_pure(:, 1), std_ga_pure(:, 1));
disp(' ');
disp('======================================================');
disp('         RESULTS TABLE: Case 1a: Series (u -> i)');
disp('======================================================');
disp(T_1a);

% --- Table 2: Case 1b ---
T_1b = create_results_table(N_axis', ...
    tiempos_clasico_pure(:, 2), std_clasico_pure(:, 2), ...
    tiempos_ga_pure(:, 2), std_ga_pure(:, 2));
disp(' ');
disp('======================================================');
disp('         RESULTS TABLE: Case 1b: Series (i -> u)');
disp('======================================================');
disp(T_1b);

% --- Table 3: Case 2a ---
T_2a = create_results_table(N_axis', ...
    tiempos_clasico_pure(:, 3), std_clasico_pure(:, 3), ...
    tiempos_ga_pure(:, 3), std_ga_pure(:, 3));
disp(' ');
disp('======================================================');
disp('         RESULTS TABLE: Case 2a: Parallel (i -> u)');
disp('======================================================');
disp(T_2a);

% --- Table 4: Case 2b ---
T_2b = create_results_table(N_axis', ...
    tiempos_clasico_pure(:, 4), std_clasico_pure(:, 4), ...
    tiempos_ga_pure(:, 4), std_ga_pure(:, 4));
disp(' ');
disp('======================================================');
disp('         RESULTS TABLE: Case 2b: Parallel (u -> i)');
disp('======================================================');
disp(T_2b);

disp(' ');
disp('======================================================');
disp('BENCHMARK COMPLETE');
disp('======================================================');


%==========================================================================
% --- HELPER FUNCTIONS (Do Not Modify) ---
%==========================================================================

% --- Table Formatting Helper ---
function T = create_results_table(N_values, mean_phasor, std_phasor, mean_ga, std_ga)
    % Create formatted strings: "Mean (Std)"
    Phasor_str = arrayfun(@(m, s) sprintf('%.3e (%.2e)', m, s), ...
        mean_phasor, std_phasor, 'UniformOutput', false);
    
    GA_Real_str = arrayfun(@(m, s) sprintf('%.3e (%.2e)', m, s), ...
        mean_ga, std_ga, 'UniformOutput', false);
    
    % Create the table
    N = N_values;
    T = table(N, Phasor_str, GA_Real_str);
    T.Properties.VariableNames = {'N', 'Phasor_Time_s', 'GA_Real_Time_s'};
end


% --- PARSING AND CONSTRUCTION FUNCTIONS (from v13.16) ---
function [coeficientes_finales] = parsear_vector_ga_ful(vector_ga_ful, dimension)
    coeficientes_finales = zeros(1, dimension);
    vector_part_texto = char(vector_ga_ful.GetVectorPart().ToString());
    vector_part_texto = strrep(vector_part_texto, ',', '.');
    if isempty(vector_part_texto) || strcmp(vector_part_texto, '0'), return; end
    vector_part_texto = strrep(vector_part_texto, ' ', '');
    vector_part_texto = strrep(vector_part_texto, '-', '+-');
    terminos = strsplit(vector_part_texto, '+');
    for i = 1:length(terminos)
        termino = terminos{i};
        if isempty(termino), continue; end
        datos = sscanf(termino, '%f<%d>');
        if length(datos) == 2
            coeficiente = datos(1);
            indice_base_0 = datos(2);
            indice_matlab_1 = indice_base_0 + 1;
            if indice_matlab_1 > 0 && indice_matlab_1 <= dimension
                coeficientes_finales(indice_matlab_1) = coeficiente;
            end
        end
    end
end
%==========================================================================
% --- NATIVE CALCULATION FUNCTIONS (Phasor and GA-Real) ---
% (These functions return 's' for a quick benchmark)
%==========================================================================
% --- Case 1a ---
function s = loop_ga_1a(n, coef_cos_u, coef_sen_u, Xk_vector, R_val, iter_bulk)
    s = 0; 
    output_vec = zeros(1, 2*n); 
    for iter = 1:iter_bulk 
        for k = 1:n 
            Ukx = coef_cos_u(k); Uky = coef_sen_u(k); Xk = Xk_vector(k);
            denominador = R_val^2 + Xk^2;
            if denominador < 1e-9
                output_vec(2*k-1) = 0; 
                output_vec(2*k)   = 0; 
            else
                output_vec(2*k-1) = (Ukx * R_val - Uky * Xk) / denominador; 
                output_vec(2*k)   = (Ukx * Xk + Uky * R_val) / denominador; 
            end
        end
        s = s + sum(output_vec); 
    end
end
function s = loop_fasor_1a(n, coef_cos_u, coef_sen_u, Xk_vector, R_val, iter_bulk)
    s = 0; 
    output_vec = zeros(1, 2*n); 
    for iter = 1:iter_bulk 
        for k = 1:n 
            Ukx = coef_cos_u(k); Uky = coef_sen_u(k); Xk = Xk_vector(k);
            Uh_phasor = Ukx - 1j * Uky;
            Zh_phasor = R_val + 1j * Xk;
            if abs(Zh_phasor) < 1e-9, Ih_phasor = 0;
            else, Ih_phasor = Uh_phasor / Zh_phasor;
            end
            output_vec(2*k-1) = real(Ih_phasor); 
            output_vec(2*k)   = -imag(Ih_phasor); 
        end
        s = s + sum(output_vec); 
    end
end
% --- Case 1b ---
function s = loop_ga_1b(n, coef_cos_i, coef_sen_i, Xk_vector, R_val, iter_bulk)
    s = 0;
    output_vec = zeros(1, 2*n); 
    for iter = 1:iter_bulk
        for k = 1:n
            Ikx = coef_cos_i(k); Iky = coef_sen_i(k); Xk = Xk_vector(k);
            output_vec(2*k-1) = (Ikx * R_val + Iky * Xk); 
            output_vec(2*k)   = (Iky * R_val - Ikx * Xk); 
        end
        s = s + sum(output_vec); 
    end
end
function s = loop_fasor_1b(n, coef_cos_i, coef_sen_i, Xk_vector, R_val, iter_bulk)
    s = 0;
    output_vec = zeros(1, 2*n); 
    for iter = 1:iter_bulk
        for k = 1:n
            Ikx = coef_cos_i(k); Iky = coef_sen_i(k); Xk = Xk_vector(k);
            Ih_phasor = Ikx - 1j * Iky;
            Zh_phasor = R_val + 1j * Xk;
            Uh_phasor = Ih_phasor * Zh_phasor;
            output_vec(2*k-1) = real(Uh_phasor); 
            output_vec(2*k)   = -imag(Uh_phasor); 
        end
        s = s + sum(output_vec); 
    end
end
% --- Case 2a ---
function s = loop_ga_2a(n, coef_cos_i, coef_sen_i, Bk_vector, G_val, iter_bulk)
    s = 0;
    output_vec = zeros(1, 2*n); 
    for iter = 1:iter_bulk
        for k = 1:n
            Ikx = coef_cos_i(k); Iky = coef_sen_i(k); Bk = Bk_vector(k);
            denominador = G_val^2 + Bk^2;
            if denominador < 1e-9
                output_vec(2*k-1) = 0; 
                output_vec(2*k)   = 0; 
            else
                output_vec(2*k-1) = (Ikx * G_val - Iky * Bk) / denominador; 
                output_vec(2*k)   = (Ikx * Bk + Iky * G_val) / denominador; 
            end
        end
        s = s + sum(output_vec); 
    end
end
function s = loop_fasor_2a(n, coef_cos_i, coef_sen_i, Bk_vector, G_val, iter_bulk)
    s = 0;
    output_vec = zeros(1, 2*n); 
    for iter = 1:iter_bulk
        for k = 1:n
            Ikx = coef_cos_i(k); Iky = coef_sen_i(k); Bk = Bk_vector(k);
            Ih_phasor = Ikx - 1j * Iky;
            Yh_phasor = G_val + 1j * Bk;
            if abs(Yh_phasor) < 1e-9, Uh_phasor = 0;
            else, Uh_phasor = Ih_phasor / Yh_phasor;
            end
            output_vec(2*k-1) = real(Uh_phasor); 
            output_vec(2*k)   = -imag(Uh_phasor); 
        end
        s = s + sum(output_vec); 
    end
end
% --- Case 2b ---
function s = loop_ga_2b(n, coef_cos_u, coef_sen_u, Bk_vector, G_val, iter_bulk)
    s = 0;
    output_vec = zeros(1, 2*n); 
    for iter = 1:iter_bulk
        for k = 1:n
            Ukx = coef_cos_u(k); Uky = coef_sen_u(k); Bk = Bk_vector(k);
            output_vec(2*k-1) = (Ukx * G_val + Uky * Bk); 
            output_vec(2*k)   = (Uky * G_val - Ukx * Bk); 
        end
        s = s + sum(output_vec); 
    end
end
function s = loop_fasor_2b(n, coef_cos_u, coef_sen_u, Bk_vector, G_val, iter_bulk)
    s = 0;
    output_vec = zeros(1, 2*n); 
    for iter = 1:iter_bulk
        for k = 1:n
            Ukx = coef_cos_u(k); Uky = coef_sen_u(k); Bk = Bk_vector(k);
            Uh_phasor = Ukx - 1j * Uky;
            Yh_phasor = G_val + 1j * Bk;
            Ih_phasor = Uh_phasor * Yh_phasor;
            output_vec(2*k-1) = real(Ih_phasor); 
            output_vec(2*k)   = -imag(Ih_phasor); 
        end
        s = s + sum(output_vec); 
    end
end