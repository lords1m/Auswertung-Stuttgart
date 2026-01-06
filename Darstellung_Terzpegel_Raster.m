%% ============================================================
%  Darstellung_Terzpegel_Raster.m
%
%  Erstellt für jede gefundene Variante einen Plot mit einem
%  4x4 Kachelfeld (Raster). In jeder Kachel wird der Terzpegel-
%  verlauf (Spektrum) der entsprechenden Position angezeigt.
%
%  Reihenfolge: Links oben nach rechts unten (Zeilenweise).
%  Layout:
%    M1  M2  M3  M4
%    M5  M6  M7  M8
%    M9  M10 M11 M12
%    Q1  M13 M14 M15
%
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

%% ---------------- Einstellungen ----------------
dataDir = 'data';
outputPlotDir = 'Plots';
fs = 500e3; % 500 kHz Abtastrate

% Normgerechte 1/3-Oktav-Mittenfrequenzen (IEC 61260)
f_terz = double([ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ]);

% Layout Definition (4x4)
positionsLayout = {
    'M1',  'M2',  'M3',  'M4';
    'M5',  'M6',  'M7',  'M8';
    'M9',  'M10', 'M11', 'M12';
    'Q1',  'M13', 'M14', 'M15';
};

%% ---------------- Setup & Dateiprüfung ----------------
if ~exist(dataDir, 'dir'), error('Datenordner "%s" nicht gefunden!', dataDir); end
if ~exist(outputPlotDir,'dir'), mkdir(outputPlotDir); end

% Alle .mat Dateien finden
dirInfo = dir(fullfile(dataDir, '*.mat'));
matFiles = {dirInfo.name};
if isempty(matFiles), error('Keine .mat-Dateien im Ordner "%s" gefunden!', dataDir); end

%% ---------------- Varianten identifizieren ----------------
% Wir suchen nach eindeutigen Variantennamen in den Dateinamen
variantNames = {};
for i = 1:numel(matFiles)
    % Versuche Muster "Variante_X..." zu erkennen
    % Ignoriere "Pos_..." Teil für den Variantennamen
    tokens = regexp(matFiles{i}, '^(.*?)[_,]Pos', 'tokens', 'once', 'ignorecase');
    if isempty(tokens)
        tokens = regexp(matFiles{i}, '^(.*?)[_,]Quelle', 'tokens', 'once', 'ignorecase');
    end
    
    if ~isempty(tokens)
        variantNames{end+1} = tokens{1};
    end
end
variantNames = unique(variantNames);
fprintf('Gefundene Varianten: %d\n', numel(variantNames));
disp(variantNames);

%% ---------------- Globale Referenz (FS) finden ----------------
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
    catch
        continue;
    end
end
if FS_global == 0
    warning('Konnte FS_global nicht bestimmen. Setze auf 1.');
    FS_global = 1;
else
    fprintf('✓ Globale Referenz (FS) ermittelt: %g\n', FS_global);
end

%% ---------------- Daten sammeln & Globales Maximum finden ----------------
fprintf('Berechne Summenpegel für alle Varianten...\n');
variantData = struct();
allValues = [];

for v = 1:numel(variantNames)
    variante = variantNames{v};
    gridData = NaN(4,4);
    
    [rows, cols] = size(positionsLayout);
    for r = 1:rows
        for c = 1:cols
            posName = positionsLayout{r, c};
            filePath = find_mat_file(dataDir, matFiles, variante, posName);
            
            if ~isempty(filePath)
                L_dBFS = calculateTerzpegelForFile(filePath, fs, f_terz, FS_global);
                if ~all(isnan(L_dBFS))
                    % Energetische Summe
                    sumLevel = 10 * log10(sum(10.^(L_dBFS / 10)));
                    gridData(r,c) = sumLevel;
                end
            end
        end
    end
    variantData(v).name = variante;
    variantData(v).grid = gridData;
    allValues = [allValues; gridData(:)];
end

% Globale Farbskala bestimmen
validVals = allValues(~isnan(allValues));
if isempty(validVals)
    cLim = [-60 0];
else
    cLim = [floor(min(validVals)), ceil(max(validVals))];
    if diff(cLim) < 5, cLim = [cLim(1)-2.5, cLim(2)+2.5]; end
end
fprintf('Farbskala gesetzt auf: %.1f bis %.1f dBFS\n', cLim(1), cLim(2));

%% ---------------- Plots erstellen ----------------
for v = 1:numel(variantNames)
    variante = variantData(v).name;
    gridData = variantData(v).grid;
    
    fprintf('Erstelle Heatmap für: %s\n', variante);
    
    fig = figure('Name', ['Summenpegel: ' variante], 'Position', [100, 100, 800, 700], 'Visible', 'off');
    
    % Heatmap
    imagesc(gridData);
    colormap(jet);
    caxis(cLim); % Absolute Skalierung
    cb = colorbar;
    cb.Label.String = 'Summenpegel [dBFS]';
    cb.Label.FontSize = 12;
    
    title(['Summenpegel - ' strrep(variante, '_', ' ')], 'FontSize', 16, 'FontWeight', 'bold');
    axis square; axis off;
    
    % Text-Labels
    [rows, cols] = size(positionsLayout);
    for r = 1:rows
        for c = 1:cols
            val = gridData(r,c);
            posName = positionsLayout{r,c};
            
            if ~isnan(val)
                str = sprintf('%s\n%.1f dB', posName, val);
                textColor = 'k'; 
            else
                str = sprintf('%s\nN/A', posName);
                textColor = [0.5 0.5 0.5];
            end
            
            text(c, r, str, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'Color', textColor);
        end
    end
    
    % Speichern
    plotDir = fullfile(outputPlotDir, variante);
    if ~exist(plotDir, 'dir'), mkdir(plotDir); end
    
    fileName = fullfile(plotDir, ['Summenpegel_Heatmap_' variante '.png']);
    saveas(fig, fileName);
    close(fig);
end

fprintf('\nFertig.\n');

%% ---------------- HILFSFUNKTIONEN ----------------

function filePath = find_mat_file(dataDir, allFiles, variante, posName)
    filePath = '';
    
    if startsWith(posName, 'M')
        posNum = extractAfter(posName, 'M');
        
        % Regex für exakten Match: Start mit Variante, dann Trenner, dann Pos+Nummer
        % (?i) = case insensitive, 0* = optionale führende Nullen
        pattern = ['^' regexptranslate('escape', variante) '(?i)[_,]Pos[_,]?0*' posNum '\.mat$'];

        for i = 1:numel(allFiles)
            fname = allFiles{i};
            if ~isempty(regexp(fname, pattern, 'once'))
                filePath = fullfile(dataDir, fname);
                return;
            end
        end
        
    elseif startsWith(posName, 'Q')
        pattern = ['^' regexptranslate('escape', variante) '(?i)[_,]Quelle\.mat$'];
        for i = 1:numel(allFiles)
            fname = allFiles{i};
            if ~isempty(regexp(fname, pattern, 'once'))
                filePath = fullfile(dataDir, fname);
                return;
            end
        end
    end
end

function L_dBFS = calculateTerzpegelForFile(filePath, fs, f_terz, FS_global)
    L_dBFS = NaN(1, length(f_terz));
    try
        S = load(filePath);
        ir = extractIR(S);
        if isempty(ir), return; end
        
        [ir_trunc, ~, ~, ~, ~, ~] = truncateIR(ir);
        N = length(ir_trunc);
        IR_fft = fft(ir_trunc, N);
        freq = (0:N-1) * (fs / N);
        
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
    catch
    end
end

function ir = extractIR(S)
    ir = [];
    if isfield(S,'RiR') && ~isempty(S.RiR), ir = double(S.RiR(:));
    elseif isfield(S,'RIR') && ~isempty(S.RIR), ir = double(S.RIR(:));
    elseif isfield(S,'aufn') && ~isempty(S.aufn), ir = double(S.aufn(:)); % Fallback
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

function [ir_trunc, start_idx, end_idx, E_ratio, SNR_dB, dynamic_range_dB] = truncateIR(ir)
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