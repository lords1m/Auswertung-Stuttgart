%% ============================================================ 
%  Terzband-Nachhallzeit (RT60) aus RIR-Dateien (Angepasst)
% 
%  Methode: Schroeder-Integration + T20
%  Struktur:
%  1. Stufe: Verarbeitet jede ausgewählte Variante, berechnet RT60 für
%     jede Position, speichert Ergebnisse in Excel und plottet
%     Einzel- und Mittelwert-RT60.
%  2. Stufe: Vergleicht die Mittelwerte aller verarbeiteten Varianten
%     in einem finalen Gesamtplot.
% 
%  Anpassungen:
%  - Lädt .mat Dateien aus zentralem /data Ordner
%  - Speichert alle Ausgaben in /Plots und /Excel
%  - Flexible Auswahl von Varianten und Positionen
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
dataDir = 'data';
outputPlotDir = 'Plots';
outputExcelDir = 'Excel';

% Varianten zum Verarbeiten (leer = automatisch alle außer '_alt' suchen)
selectedVariants = {'Variante_4'};  % z.B. {'Variante_1', 'Variante_2'} oder {} für Auto
excludePattern = '_alt'; % Ausschlussmuster für Varianten

% Positionen (Messpunkte) auswählen
selectedPositions = 1:14;   % z.B. [1 3 5] oder 1:14
fs = 500e3;                 % Abtastrate [Hz]

% Positionen, die für den Mittelwert pro Variante verwendet werden
positions_to_average = [5, 6, 7, 9, 10, 11, 13, 14];

% Normgerechte 1/3-Oktav-Mittenfrequenzen (IEC 61260)
f_terz = double([ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ]);
nFreq = numel(f_terz);

%% ---------------- Setup & Varianten ermitteln ----------------
if ~exist(dataDir, 'dir'), error('Datenordner "%s" nicht gefunden!', dataDir); end
if ~exist(outputPlotDir,'dir'), mkdir(outputPlotDir); end
if ~exist(outputExcelDir,'dir'), mkdir(outputExcelDir); end

dirInfo = dir(fullfile(dataDir, '*.mat'));
matFiles = {dirInfo.name};
if isempty(matFiles), error('Keine .mat-Dateien im Ordner "%s" gefunden!', dataDir); end

fprintf('Gefundene .mat-Dateien: %d\n', numel(matFiles));

% Varianten aus Dateinamen extrahieren
variantNamesRaw = cell(numel(matFiles), 1);
validIdx = false(numel(matFiles), 1);
for i = 1:numel(matFiles)
    tokens = regexp(matFiles{i}, '^([^,]+?)(?:_neu)?[,_]Pos', 'tokens', 'once');
    if ~isempty(tokens)
        variantName = tokens{1};
        if isempty(excludePattern) || ~contains(variantName, excludePattern)
            variantNamesRaw{i} = variantName;
            validIdx(i) = true;
        end
    end
end
variantNames = unique(variantNamesRaw(validIdx));

if ~isempty(selectedVariants)
    variantNames = intersect(variantNames, selectedVariants, 'stable');
end

if isempty(variantNames), error('Keine Varianten zum Verarbeiten gefunden (nach Filter).'); end

fprintf('\nVerarbeite Varianten:\n');
fprintf('  - %s\n', variantNames{:});

% Erstelle File-Lookup-Map für schnelleren Zugriff
fileMap = buildFileMap(matFiles, variantNames, selectedPositions);


%% ========================================================================
%  Stufe 1: Verarbeitung pro Variante
% ========================================================================

for vi = 1:numel(variantNames)
    variantName = variantNames{vi};
    fprintf('\n========================================\n');
    fprintf('Verarbeite %s...\n', variantName);
    fprintf('========================================\n');

    nPos = numel(selectedPositions);
    RT60_data = NaN(nPos, nFreq);   % [Messung x Terzband] 

    for pi = 1:nPos
        pos = selectedPositions(pi);
        fname = getFileFromMap(fileMap, variantName, pos);
        
        if isempty(fname)
            warning('Keine Datei für Pos %d in %s gefunden. Position übersprungen.', pos, variantName);
            continue;
        end

        fprintf('  Verarbeite: %s, Pos %d\n', variantName, pos);
        
        S = load(fullfile(dataDir, fname));
        ir = extractIR(S);
        
        if isempty(ir)
            warning('Datei %s enthält keinen erkennbaren IR-Vektor. Übersprungen.', fname);
            continue;
        end

        % RT60-Berechnung für jedes Terzband
        for k = 1:nFreq
            RT60_data(pi,k) = calculate_rt60_band(ir, f_terz(k), fs);
        end
    end

    % --- Excel-Export ---
    rowNames = compose('Pos_%02d', selectedPositions);
    colNames = compose('F%.0f', f_terz);
    T_RT60 = array2table(RT60_data, 'RowNames', rowNames, 'VariableNames', colNames);
    excelFile = fullfile(outputExcelDir, sprintf('RT60_%s.xlsx', variantName));
    try
        writetable(T_RT60, excelFile, 'WriteRowNames', true);
        fprintf('Excel geschrieben: %s\n', excelFile);
    catch ME
        warning(ME.identifier, 'Konnte Excel nicht schreiben: %s', ME.message);
    end

    % --- Plots pro Variante ---
    y_min_RT60 = min(RT60_data(:), [], 'omitnan');
    y_max_RT60 = max(RT60_data(:), [], 'omitnan');
    y_range_RT60 = [floor(y_min_RT60*10)/10, ceil(y_max_RT60*10)/10];
    if isempty(y_range_RT60) || any(isinf(y_range_RT60)) || any(isnan(y_range_RT60)), y_range_RT60 = [0 1]; end
    
    xtick_vals = [4000 5000 10000 20000 50000 60000];
    xtick_labels = {'4k', '5k', '10k', '20k', '50k', '60k'};
    f_lim = [4000 60000];
    
    % Plot 1: RT60 pro Position
    for pi = 1:nPos
        pos = selectedPositions(pi);
        if all(isnan(RT60_data(pi, :))), continue; end
        
        fig = figure('Visible','off', 'Position', [100, 100, 1000, 500]);
        stairs(f_terz, RT60_data(pi,:), 'LineWidth', 2);
        grid on;
        set(gca, 'XScale', 'log');
        xlabel('Frequenz [Hz]'); ylabel('Nachhallzeit RT60 [s]');
        title(sprintf('RT60 - %s - Position %02d', variantName, pos));
        xlim(f_lim); ylim(y_range_RT60);
        xticks(xtick_vals); xticklabels(xtick_labels);

        filename = fullfile(outputPlotDir, sprintf('RT60_%s_Pos_%02d', variantName, pos));
        saveas(fig, [filename '.png']);
        close(fig);
    end

    % Plot 2: Mittelwert ausgewählter Positionen
    % Stellen sicher, dass die zu mittelnden Positionen auch in selectedPositions enthalten sind
    avg_pos_indices = find(ismember(selectedPositions, positions_to_average));
    if ~isempty(avg_pos_indices)
        RT60_mean = mean(RT60_data(avg_pos_indices, :), 1, 'omitnan');

        fig_mean = figure('Visible','off', 'Position', [100, 100, 1000, 500]);
        stairs(f_terz, RT60_mean, 'LineWidth', 2.5, 'Color', [0.8500 0.3250 0.0980]);
        grid on;
        set(gca, 'XScale', 'log');
        xlabel('Frequenz [Hz]'); ylabel('Nachhallzeit RT60 [s]');
        title(sprintf('%s - Gemittelte RT60', strrep(variantName,'_',' ')));
        xlim(f_lim); ylim(y_range_RT60);
        xticks(xtick_vals); xticklabels(xtick_labels);
        
        filename_mean = fullfile(outputPlotDir, sprintf('RT60_%s_Mittelwert', variantName));
        saveas(fig_mean, [filename_mean '.png']);
        close(fig_mean);
        fprintf('Mittelwert-Plot erstellt: %s.png\n', filename_mean);
    end
    fprintf('%s abgeschlossen!\n', variantName);
end

%% ========================================================================
%  Stufe 2: Vergleich der Mittelwerte aller verarbeiteten Varianten
% ========================================================================
fprintf('\n========================================\n');
fprintf('Lade RT60-Daten für Variantenvergleich...\n');
fprintf('========================================\n');

RT60_means = NaN(numel(variantNames), nFreq);

for vi = 1:numel(variantNames)
    varName = variantNames{vi};
    excelPath = fullfile(outputExcelDir, sprintf('RT60_%s.xlsx', varName));

    if exist(excelPath, 'file')
        T = readtable(excelPath, 'ReadRowNames', true);
        
        % Finde die Indizes der zu mittelnden Positionen in der Tabelle
        posInTable = T.Properties.RowNames; % z.B. {'Pos_01', ...}
        posToAvgStr = compose('Pos_%02d', positions_to_average);
        [~, avg_indices] = intersect(posInTable, posToAvgStr);
        
        if ~isempty(avg_indices)
            RT60_matrix = table2array(T);
            RT60_means(vi, :) = mean(RT60_matrix(avg_indices, :), 1, 'omitnan');
            fprintf('✓ %s geladen und gemittelt\n', varName);
        else
            warning('Für %s konnten keine der Positionen für den Mittelwert in der Excel-Datei gefunden werden.', varName);
        end
    else
        warning('Excel-Datei nicht gefunden für den Vergleich: %s', excelPath);
    end
end

% Plot erstellen: Vergleich aller Varianten
fig_compare = figure('Visible','on', 'Position', [100, 100, 1200, 600]);
hold on;
for vi = 1:numel(variantNames)
    if ~all(isnan(RT60_means(vi,:)))
        stairs(f_terz, RT60_means(vi, :), 'LineWidth', 2.5, 'DisplayName', strrep(variantNames{vi}, '_', ' '));
    end
end
hold off;

grid on;
set(gca, 'XScale', 'log');
xlabel('Frequenz [Hz]', 'FontSize', 12);
ylabel('Nachhallzeit RT60 [s]', 'FontSize', 12);
title('Vergleich der RT60-Mittelwerte aller Varianten', 'FontSize', 14);
xlim(f_lim);

 y_min_all = min(RT60_means(:), [], 'omitnan');
 y_max_all = max(RT60_means(:), [], 'omitnan');
 if ~isempty(y_min_all) && ~isempty(y_max_all), ylim([floor(y_min_all*10)/10, ceil(y_max_all*10)/10]); end

xticks(xtick_vals); xticklabels(xtick_labels);
legend('Location', 'best', 'FontSize', 11);

filename_compare = fullfile(outputPlotDir, 'Vergleich_RT60_Mittelwerte_Alle_Varianten');
saveas(fig_compare, [filename_compare '.png']);
saveas(fig_compare, [filename_compare '.fig']);

fprintf('\n========================================\n');
fprintf('Vergleichsplot erstellt: %s.png\n', filename_compare);
fprintf('========================================\n');

close(fig_compare);
disp('Alle Varianten erfolgreich verarbeitet!');

%% ========================================================================
%  HILFSFUNKTIONEN
% ========================================================================

function rt60 = calculate_rt60_band(ir, f_center, fs)
    filterOrder = 4;
    rt60 = NaN;

    % Terzband-Filter
    f1 = f_center / 2^(1/6);
    f2 = f_center * 2^(1/6);
    Wn = [f1 f2] / (fs/2);
    if Wn(2) >= 1, return; end
    
    % Umstellung auf SOS (Second-Order Sections) zur Vermeidung numerischer Instabilität
    [sos, g] = butter(filterOrder/2, Wn, 'bandpass');
    
    % Einfache Truncation für die spezifische Berechnung
    [ir_trunc, ~, ~, ~, ~, ~] = truncateIR(ir, fs);
    
    % Filterung mit SOS-Format
    ir_filt = filtfilt(sos, g, ir_trunc);
    t = (0:length(ir_trunc)-1)'/fs;

    % Schroeder-Integration
    E = flipud(cumsum(flipud(ir_filt.^2)));
    if max(E) == 0, return; end
    E_dB = 10*log10(E / max(E) + eps);

    % T20-Bereich für die Regression
    idx_upper = find(E_dB <= -5, 1, 'first');
    idx_lower = find(E_dB <= -25, 1, 'first');
    if isempty(idx_upper) || isempty(idx_lower) || (idx_lower - idx_upper < 10), return; end
    idx = idx_upper:idx_lower;

    % Lineare Regression
    p = polyfit(t(idx), E_dB(idx), 1);
    if p(1) >= 0, return; end % Steigung muss negativ sein

    rt60 = -60 / p(1);
end


function fileMap = buildFileMap(matFiles, variantNames, positions)
    nVariants = numel(variantNames);
    maxPos = max(positions);
    fileMap = cell(nVariants, maxPos);
    
    for i = 1:numel(matFiles)
        fname = matFiles{i};
        
        for vi = 1:nVariants
            if startsWith(fname, variantNames{vi})
                tok = regexp(fname, 'Pos[_,]?(\d+)', 'tokens', 'once');
                if ~isempty(tok)
                    posNum = str2double(tok{1});
                    if ismember(posNum, positions)
                        fileMap{vi, posNum} = fname;
                    end
                end
                break;
            end
        end
    end
end

function fname = getFileFromMap(fileMap, variantName, pos)
    [nVariants, ~] = size(fileMap);
    
    for vi = 1:nVariants
       if ~isempty(fileMap{vi,pos})
           % Quick check based on variant name start - assumes buildFileMap order
           if startsWith(fileMap{vi,pos}, variantName)
               fname = fileMap{vi,pos};
               return;
           end
       end
    end
    fname = ''; % Fallback if not found
end

function ir = extractIR(S)
    ir = [];
    if isfield(S,'RiR') && ~isempty(S.RiR), ir = double(S.RiR(:));
    elseif isfield(S,'RIR') && ~isempty(S.RIR), ir = double(S.RIR(:));
    else
        fns = fieldnames(S);
        for f = 1:numel(fns)
            fname = fns{f};
            if startsWith(fname, '__'), continue; end
            v = S.(fname);
            if isnumeric(v) && numel(v) > 1000, ir = double(v(:)); return; end
        end
    end
    if ~isempty(ir) && numel(ir) < 2, ir = []; end
end

function [ir_trunc, start_idx, end_idx, E_ratio, SNR_dB, dynamic_range_dB] = truncateIR(ir, fs)
    N_original = length(ir);
    ir_abs = abs(ir);
    max_amp = max(ir_abs);
    if max_amp == 0, ir_trunc = ir; start_idx=1; end_idx=N_original; E_ratio=100; SNR_dB=0; dynamic_range_dB=0; return; end

    noise_samples = ir_abs(end-round(N_original*0.1):end);
    noise_level = mean(noise_samples);
    noise_rms = std(noise_samples);

    threshold = max((noise_level + 3*noise_rms) * 10, max_amp * 0.001);
    sig_idx = find(ir_abs > threshold, 1, 'last');

    if isempty(sig_idx) || sig_idx < 100, end_idx = N_original;
    else, end_idx = min(N_original, round(sig_idx * 1.2)); end

    start_threshold = max_amp * 0.05;
    start_idx = find(ir_abs > start_threshold, 1, 'first');
    if isempty(start_idx) || start_idx < 1, start_idx = 1; end

    ir_trunc = ir(start_idx:end_idx);
    
    E_original = sum(ir.^2);
    E_truncated = sum(ir_trunc.^2);
    if E_original == 0, E_ratio = 100; else, E_ratio = E_truncated / E_original * 100; end

    signal_rms = sqrt(mean(ir_trunc.^2));
    SNR_dB = 20*log10(signal_rms / (noise_rms + eps));
    dynamic_range_dB = 20*log10(max_amp / (noise_level + eps));
end