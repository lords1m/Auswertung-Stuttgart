%% step2_compare_variants.m
%  Vergleich zweier Varianten basierend auf Processed Data
clear; clc; close all;
addpath('functions');

% --- Einstellungen ---
procDir = 'processed';
var1_name = 'Variante_1';
var2_name = 'Variante_2';
pos_to_compare = 13; 

% Laden der Zusammenfassung (optional, um Dateinamen zu finden)
load(fullfile(procDir, 'Summary_Database.mat'), 'summary_table');

% --- Daten laden (über Helper oder direkt via Dateinamen-Konvention) ---
file1 = fullfile(procDir, sprintf('Proc_%s_Pos%d.mat', var1_name, pos_to_compare));
file2 = fullfile(procDir, sprintf('Proc_%s_Pos%d.mat', var2_name, pos_to_compare));

if ~exist(file1, 'file') || ~exist(file2, 'file')
    error('Dateien nicht gefunden. Bitte step1_process_data.m ausführen.');
end

D1 = load(file1); R1 = D1.Result;
D2 = load(file2); R2 = D2.Result;

% --- Plot: Terzpegel ---
f = R1.freq.f_center;

figure('Color', 'w', 'Position', [100 100 1000 500]);

% Spektren
subplot(1,2,1);
semilogx(f, R1.freq.terz_dbfs, 'b-o', 'LineWidth', 1.5, 'DisplayName', var1_name);
hold on;
semilogx(f, R2.freq.terz_dbfs, 'r-x', 'LineWidth', 1.5, 'DisplayName', var2_name);
grid on;
xlabel('Frequenz [Hz]'); ylabel('Pegel [dBFS]');
title(sprintf('Vergleich Position %d', pos_to_compare));
legend show;
xlim([4000 60000]);

% Differenz
subplot(1,2,2);
diff_dB = R2.freq.terz_dbfs - R1.freq.terz_dbfs;
bar(categorical(f), diff_dB); % Balkendiagramm oft besser für Terzen
title('Differenz (Var2 - Var1)');
ylabel('Delta [dB]');
grid on;

% --- Ausgabe Metriken ---
fprintf('Vergleich Position %d:\n', pos_to_compare);
fprintf('%s Summenpegel: %.2f dBFS (SNR: %.1f dB)\n', var1_name, R1.freq.sum_level, R1.time.metrics.snr_db);
fprintf('%s Summenpegel: %.2f dBFS (SNR: %.1f dB)\n', var2_name, R2.freq.sum_level, R2.time.metrics.snr_db);