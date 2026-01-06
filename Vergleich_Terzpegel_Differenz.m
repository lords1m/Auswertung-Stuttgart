%% ============================================================
%  Vergleich von zwei Terzpegel-Messungen (oder Gruppen)
%  Ergebnis: Plot der Differenz in dB
%  Update: Kann nun den Mittelwert aus mehreren Varianten/Positionen
%          bilden und diese Gruppen dann vergleichen.
% ============================================================

clear;
clc;
close all;

% Stelle sicher, dass wir im richtigen Verzeichnis sind
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir)
    cd(scriptDir);
end
fprintf('Arbeitsverzeichnis: %s\n', pwd);

%% ---------------- Einstellungen (anpassen) ----------------
% --- Messung 1 (Referenz) ---
% Kann eine oder mehrere Varianten/Positionen sein.
% Bei mehreren wird der energetische Mittelwert gebildet.
variante1 = {'Variante_1', 'Variante_3', 'Variante_4'}; % z.B. {'Variante_1', 'Variante_2'}
position1 = 7:14;              % z.B. 1:14 oder [1, 3, 5]

% --- Messung 2 (Vergleich) ---
variante2 = {'Variante_2'};
position2 = 7:14;

% --- Allgemeine Einstellungen ---
dataDir = 'data'; % Datenordner mit allen .mat-Dateien
fs = 500e3; % 500 kHz Abtastrate

% Normgerechte 1/3-Oktav-Mitttenfrequenzen (IEC 61260)
f_terz = double([ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ]);

%% ---------------- Prüfe data-Ordner ----------------
if ~exist(dataDir, 'dir')
    error('Datenordner "%s" nicht gefunden!', dataDir);
end

%% ---------------- Dateien einlesen und globale Referenz finden ----------------
dirInfo = dir(fullfile(dataDir, '*.mat'));
matFiles = {dirInfo.name};
if isempty(matFiles)
    error('Keine .mat-Dateien im Ordner "%s" gefunden!', dataDir);
end

FS_global = 0;
fprintf('Suche globale Referenz (FS) in allen .mat Dateien...\n');
for i = 1:numel(matFiles)
    try
        S = load(fullfile(dataDir, matFiles{i}));
        ir = extractIR(S);
        if ~isempty(ir)
            [ir_trunc, ~, ~, ~, ~, ~] = truncateIR(ir);
            N = numel(ir_trunc);
            H_mag = abs(fft(ir_trunc, N));
            FS_global = max(FS_global, max(H_mag(1:floor(N/2)+1)));
        end
    catch ME
        warning('Datei %s konnte nicht für FS-Bestimmung gelesen werden: %s', matFiles{i}, ME.message);
    end
end

if FS_global == 0
    error('Keine gültigen Impulsantworten gefunden. FS_Global konnte nicht bestimmt werden.');
end
fprintf('✓ Globale Referenz (FS) ermittelt: %g\n', FS_global);


%% ---------------- Terzpegel für beide Messungen/Gruppen berechnen ----------------
fprintf('\nBerechne Terzpegel für Gruppe 1...\n');
terzpegel1 = calculate_mean_terzpegel(variante1, position1, dataDir, matFiles, fs, f_terz, FS_global);

fprintf('\nBerechne Terzpegel für Gruppe 2...\n');
terzpegel2 = calculate_mean_terzpegel(variante2, position2, dataDir, matFiles, fs, f_terz, FS_global);

if all(isnan(terzpegel1)) || all(isnan(terzpegel2))
    error('Mindestens eine der Messgruppen konnte nicht ausgewertet werden. Abbruch.');
end

%% ---------------- Differenz berechnen und plotten ----------------
differenz_dB = terzpegel2 - terzpegel1;

% Beschreibungen für den Plot erstellen
desc1 = create_description_string(variante1, position1);
desc2 = create_description_string(variante2, position2);

fig = figure('Position',[100, 100, 1200, 600]);
stairs(1:length(f_terz), differenz_dB, 'b-', 'LineWidth', 2);
hold on;
yline(0, 'k--', 'LineWidth', 1.5); % Nulllinie für Referenz
grid on;

% Achsen und Titel
title(sprintf('Differenz Terzpegel: (%s) - (%s)', desc2, desc1), 'FontSize', 14);
xlabel('Frequenz [Hz]', 'FontSize', 12);
ylabel('Pegeldifferenz [dB]', 'FontSize', 12);

% Ticks für die Frequenzachse setzen
set(gca, 'XTick', 1:2:length(f_terz));
set(gca, 'XTickLabel', f_terz(1:2:length(f_terz)));
xtickangle(45);
set(gca, 'XScale', 'log'); % Logarithmische Skala für Frequenz

xlim([find(f_terz > 4000, 1, 'first')-1, find(f_terz < 60000, 1, 'last')+1]);

% Y-Achsen-Limits sinnvoll setzen
max_abs_diff = max(abs(differenz_dB), [], 'omitnan');
if isempty(max_abs_diff) || max_abs_diff == 0
    ylim([-5 5]);
else
    ylim([-ceil(max_abs_diff/5)*5, ceil(max_abs_diff/5)*5]);
end

legend(sprintf('Differenz: (%s) - (%s)', desc2, desc1), 'Location', 'northwest');
set(gcf, 'Color', 'w'); % Weißer Hintergrund

fprintf('\nPlot der Differenz wurde erstellt.\n');

%% ---------------- Hilfsfunktion zur Mittelwertbildung ----------------
function mean_L_dBFS = calculate_mean_terzpegel(varianten, positionen, dataDir, allFiles, fs, f_terz, FS_global)
    
    if ischar(varianten)
        varianten = {varianten}; % Sicherstellen, dass es immer ein Cell-Array ist
    end
    
    all_L_dBFS = []; % Matrix zum Sammeln aller Pegel-Vektoren
    
    for i = 1:numel(varianten)
        variante = varianten{i};
        for j = 1:numel(positionen)
            position = positionen(j);
            
            fprintf('  Lese ein: %s, Pos %d\n', variante, position);
            
            L_dBFS_single = getTerzpegel(variante, position, dataDir, allFiles, fs, f_terz, FS_global);
            
            if ~all(isnan(L_dBFS_single))
                all_L_dBFS = [all_L_dBFS; L_dBFS_single];
            else
                warning('Keine gültigen Daten für %s, Pos %d erhalten.', variante, position);
            end
        end
    end
    
    if isempty(all_L_dBFS)
        warning('Keine einzige gültige Messung in der Gruppe gefunden.');
        mean_L_dBFS = NaN(1, length(f_terz));
        return;
    end
    
    % Energetische Mittelung
    fprintf('Bilde energetischen Mittelwert aus %d Messung(en).\n', size(all_L_dBFS, 1));
    linear_values = 10.^(all_L_dBFS / 10);
    mean_linear = mean(linear_values, 1);
    mean_L_dBFS = 10 * log10(mean_linear);
end


%% ---------------- Hilfsfunktion zur Terzpegelberechnung (Einzelmessung) ----------------
function L_dBFS = getTerzpegel(variante, position, dataDir, allFiles, fs, f_terz, FS_global)
    
    L_dBFS = NaN(1, length(f_terz)); % Initialisieren mit NaN

    % Passende Datei finden
    posStr = num2str(position);
    fname = '';
    
    % Suche nach exaktem Match zuerst (effizienter)
    searchPattern_exact = sprintf('%s_Pos%s.mat', variante, posStr);
    for i = 1:numel(allFiles)
        if strcmp(allFiles{i}, searchPattern_exact)
            fname = allFiles{i};
            break;
        end
    end
    
    % Wenn nicht gefunden, flexiblere Suche mit RegEx
    if isempty(fname)
        pattern = ['^' variante '.*Pos[_,]?' posStr '([^\d]|\.mat$)'];
        for i = 1:numel(allFiles)
             if ~isempty(regexp(allFiles{i}, pattern, 'once'))
                fname = allFiles{i};
                break;
            end
        end
    end

    if isempty(fname)
        % Keine Warnung hier, wird in der aufrufenden Funktion behandelt
        return;
    end
    
    % Lade und verarbeite Datei
    try
        S = load(fullfile(dataDir, fname));
        ir = extractIR(S);
        
        if isempty(ir)
            warning('Datei %s enthält keinen erkennbaren IR-Vektor.', fname);
            return;
        end
        
        % Lundeby-Truncation
        [ir_trunc, ~, ~, ~, ~, ~] = truncateIR(ir);
        N = length(ir_trunc);

        % FFT berechnen
        IR_fft = fft(ir_trunc, N);
        freq = (0:N-1) * (fs / N);
        
        % Terzband-Energie berechnen
        f_lower = f_terz / 2^(1/6);
        f_upper = f_terz * 2^(1/6);
        IR_fft_abs2 = abs(IR_fft).^2;

        for k = 1:length(f_terz)
            if f_upper(k) >= fs/2, continue; end
            idx = (freq >= f_lower(k)) & (freq <= f_upper(k));
            if ~any(idx), continue; end

            total_energy = sum(IR_fft_abs2(idx)) / N;
            rms_band = sqrt(total_energy);
            L_dBFS(k) = 20*log10((rms_band + eps) / FS_global);
        end
        
    catch ME
        warning('Fehler bei der Verarbeitung von %s: %s', fname, ME.message);
        return;
    end
end

%% ---------------- Plot-Beschriftung ----------------
function desc = create_description_string(varianten, positionen)
    if ischar(varianten)
        varianten = {varianten};
    end
    
    var_str = strjoin(varianten, '/');
    
    if numel(positionen) > 3
        pos_str = sprintf('Pos %d-%d', min(positionen), max(positionen));
    else
        pos_str = ['Pos ' strjoin(arrayfun(@num2str, positionen, 'UniformOutput', false), ',')];
    end
    
    if numel(varianten) > 1 || numel(positionen) > 1
        desc = sprintf('Mittelw. %s, %s', var_str, pos_str);
    else
        desc = sprintf('%s, %s', var_str, pos_str);
    end
end


%% ---------------- Hilfsfunktionen aus Terzpegel_DBFs.m ----------------
function ir = extractIR(S)
    ir = [];
    if isfield(S,'RiR') && ~isempty(S.RiR)
        ir = double(S.RiR(:));
    elseif isfield(S,'RIR') && ~isempty(S.RIR)
        ir = double(S.RIR(:));
    else
        fns = fieldnames(S);
        for f = 1:numel(fns)
            fname = fns{f};
            if startsWith(fname, '__'), continue; end
            v = S.(fname);
            if isnumeric(v) && numel(v) > 1000
                ir = double(v(:));
                return;
            end
        end
    end
    if ~isempty(ir) && numel(ir) < 2
        ir = [];
    end
end

function [ir_trunc, start_idx, end_idx, E_ratio, SNR_dB, dynamic_range_dB] = truncateIR(ir)
    N_original = length(ir);
    ir_abs = abs(ir);
    max_amp = max(ir_abs);

    noise_samples = ir_abs(end-round(N_original*0.1):end);
    noise_level = mean(noise_samples);
    noise_rms = std(noise_samples);

    threshold = max((noise_level + 3*noise_rms) * 10, max_amp * 0.001);
    sig_idx = find(ir_abs > threshold, 1, 'last');

    if isempty(sig_idx) || sig_idx < 100
        end_idx = N_original;
    else
        end_idx = min(N_original, round(sig_idx * 1.2));
    end

    start_threshold = max_amp * 0.05;
    start_idx = find(ir_abs > start_threshold, 1, 'first');
    if isempty(start_idx) || start_idx < 1
        start_idx = 1;
    end

    ir_trunc = ir(start_idx:end_idx);
    
    E_original = sum(ir.^2);
    E_truncated = sum(ir_trunc.^2);
    E_ratio = E_truncated / E_original * 100;

    signal_rms = sqrt(mean(ir_trunc.^2));
    SNR_dB = 20*log10(signal_rms / (noise_rms + eps));

    dynamic_range_dB = 20*log10(max_amp / (noise_level + eps));
end
