%% ============================================================
%  Darstellung_Heatmap_Video.m
%
%  Erstellt zeitabhängige Heatmap-Videos aus den Impulsantworten.
%  Zeigt die Energieverteilung (Summenpegel) im 4x4-Raster über
%  die Zeit (Abklingvorgang).
%
%  Output: .mp4 Dateien im Ordner 'Videos'
% ============================================================

clear;
clc;
close all;

% Arbeitsverzeichnis
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir), cd(scriptDir); end
fprintf('Arbeitsverzeichnis: %s\n', pwd);

%% ---------------- Einstellungen ----------------
dataDir = 'data';
outputVideoDir = 'Videos';
fs = 500e3; % 500 kHz

% Video-Parameter
videoFPS = 20;          % Bilder pro Sekunde im Video
timeStep_ms = 0.01;        % Zeitschritt zwischen Frames in der IR (5 ms für feinere Auflösung)
windowSize_ms = 0.01;     % Integrationsfenster für RMS (Glättung)
maxDuration_s = 0.1;    % Maximale Dauer der Analyse pro IR (um Rauschen am Ende zu ignorieren)
cLim = [-60 0];         % Farbskala in dBFS (festgelegt für Vergleichbarkeit)

% Layout (4x4)
positionsLayout = {
    'M1',  'M2',  'M3',  'M4';
    'M5',  'M6',  'M7',  'M8';
    'M9',  'M10', 'M11', 'M12';
    'Q1',  'M13', 'M14', 'M15';
};

%% ---------------- Setup ----------------
if ~exist(dataDir, 'dir'), error('Datenordner "%s" nicht gefunden!', dataDir); end
if ~exist(outputVideoDir,'dir'), mkdir(outputVideoDir); end

dirInfo = dir(fullfile(dataDir, '*.mat'));
matFiles = {dirInfo.name};
if isempty(matFiles), error('Keine .mat-Dateien gefunden!'); end

% Varianten finden
variantNames = {};
for i = 1:numel(matFiles)
    tokens = regexp(matFiles{i}, '^(.*?)[_,]Pos', 'tokens', 'once', 'ignorecase');
    if isempty(tokens), tokens = regexp(matFiles{i}, '^(.*?)[_,]Quelle', 'tokens', 'once', 'ignorecase'); end
    if ~isempty(tokens), variantNames{end+1} = tokens{1}; end
end
variantNames = unique(variantNames);
fprintf('Gefundene Varianten: %d\n', numel(variantNames));

%% ---------------- Globale Referenz (Max Amplitude) ----------------
MaxAmp_global = 0;
fprintf('Ermittle globale Referenz (Max Amplitude)...\n');
for i = 1:numel(matFiles)
    try
        S = load(fullfile(dataDir, matFiles{i}));
        ir = extractIR(S);
        if ~isempty(ir), MaxAmp_global = max(MaxAmp_global, max(abs(ir))); end
    catch, end
end
if MaxAmp_global == 0, MaxAmp_global = 1; end
fprintf('MaxAmp_global: %g\n', MaxAmp_global);

%% ---------------- Video-Erstellung ----------------
for v = 1:numel(variantNames)
    variante = variantNames{v};
    fprintf('\nVerarbeite Variante %d/%d: %s\n', v, numel(variantNames), variante);
    
    % 1. Alle IRs für diese Variante laden und vorbereiten
    [rows, cols] = size(positionsLayout);
    irs = cell(rows, cols);
    
    for r = 1:rows
        for c = 1:cols
            posName = positionsLayout{r, c};
            filePath = find_mat_file(dataDir, matFiles, variante, posName);
            
            if ~isempty(filePath)
                try
                    S = load(filePath);
                    rawIR = extractIR(S);
                    % Truncate richtet den Start (Direktschall) auf Index 1 aus
                    [ir_trunc, ~, ~, ~, ~, ~] = truncateIR(rawIR);
                    irs{r,c} = ir_trunc;
                catch
                    irs{r,c} = [];
                end
            else
                irs{r,c} = [];
            end
        end
    end
    
    % 2. Zeitvektor definieren
    % Wir nehmen an, dass alle IRs bei t=0 (Direktschall) beginnen dank truncateIR
    dt = timeStep_ms / 1000;
    winLen = windowSize_ms / 1000;
    winSamples = round(winLen * fs);
    
    % Bestimme maximale Länge für das Video
    currentMaxLen = 0;
    for r=1:rows, for c=1:cols, currentMaxLen = max(currentMaxLen, length(irs{r,c})); end, end
    duration = min(maxDuration_s, currentMaxLen / fs);
    timePoints = 0:dt:duration;
    
    if isempty(timePoints)
        warning('Keine gültigen Daten für %s. Überspringe.', variante);
        continue;
    end
    
    % 3. Video initialisieren
    videoName = fullfile(outputVideoDir, ['Heatmap_' variante '.mp4']);
    vObj = VideoWriter(videoName, 'MPEG-4');
    vObj.FrameRate = videoFPS;
    open(vObj);
    
    fig = figure('Visible', 'off', 'Position', [100, 100, 600, 500]);
    
    % Plot einmalig initialisieren (Objekte erstellen)
    hImg = imagesc(NaN(rows, cols));
    colormap(jet);
    caxis(cLim);
    colorbar;
    axis square; axis off;
    hTitle = title('', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Text-Objekte vorerstellen (leere Platzhalter)
    hText = gobjects(rows, cols);
    for r=1:rows, for c=1:cols, hText(r,c) = text(c, r, '', 'HorizontalAlignment', 'center', 'FontSize', 10); end, end
    
    nFrames = length(timePoints);
    fprintf('  Erstelle Frames (%d Schritte)...\n', nFrames);
    reverseStr = '';
    tStart = tic;
    
    % 4. Frames generieren
    for t_idx = 1:nFrames
        t_start = timePoints(t_idx);
        t_end = t_start + winLen;
        
        idx_start = round(t_start * fs) + 1;
        idx_end = round(t_end * fs);
        
        gridData = NaN(rows, cols);
        
        for r = 1:rows
            for c = 1:cols
                ir = irs{r,c};
                if isempty(ir) || idx_start > length(ir)
                    val = NaN;
                else
                    curr_idx_end = min(length(ir), idx_end);
                    segment = ir(idx_start:curr_idx_end);
                    
                    % RMS im Fenster berechnen
                    rms_val = sqrt(mean(segment.^2));
                    val = 20 * log10((rms_val + eps) / MaxAmp_global);
                end
                gridData(r,c) = val;
            end
        end
        
        % Vorhandene Objekte aktualisieren (statt neu erstellen) -> Viel schneller!
        set(hImg, 'CData', gridData);
        set(hTitle, 'String', sprintf('%s\nZeit: %.3f s', strrep(variante,'_',' '), t_start));
          
        % Text Labels
        for r = 1:rows
            for c = 1:cols
                val = gridData(r,c);
                pos = positionsLayout{r,c};
                if ~isnan(val) && val > cLim(1)
                    set(hText(r,c), 'String', sprintf('%s\n%.0f', pos, val), 'Color', 'k', 'Visible', 'on');
                else
                    set(hText(r,c), 'Visible', 'off');
                end
            end
        end
        
        frame = getframe(fig);
        writeVideo(vObj, frame);
        
        % Fortschrittsanzeige im Command Window
        if mod(t_idx, 5) == 0 || t_idx == nFrames
            percent = t_idx / nFrames * 100;
            
            elapsed = toc(tStart);
            remTime = (elapsed / t_idx) * (nFrames - t_idx);
            
            msg = sprintf('    Fortschritt: %3.0f%% (%d/%d) - Restzeit: %02.0f:%02.0f', ...
                          percent, t_idx, nFrames, floor(remTime/60), mod(remTime, 60));
            fprintf([reverseStr, msg]);
            reverseStr = repmat('\b', 1, length(msg));
        end
    end
    fprintf('\n');
    
    close(vObj);
    close(fig);
    fprintf('  Video gespeichert: %s\n', videoName);
end
fprintf('\nFertig.\n');

%% ---------------- HILFSFUNKTIONEN ----------------
function filePath = find_mat_file(dataDir, allFiles, variante, posName)
    filePath = '';
    if startsWith(posName, 'M')
        posNum = extractAfter(posName, 'M');
        pattern = ['^' regexptranslate('escape', variante) '(?i)[_,]Pos[_,]?0*' posNum '\.mat$'];
        idx = find(~cellfun(@isempty, regexp(allFiles, pattern, 'once')), 1);
        if ~isempty(idx), filePath = fullfile(dataDir, allFiles{idx}); end
    elseif startsWith(posName, 'Q')
        pattern = ['^' regexptranslate('escape', variante) '(?i)[_,]Quelle\.mat$'];
        idx = find(~cellfun(@isempty, regexp(allFiles, pattern, 'once')), 1);
        if ~isempty(idx), filePath = fullfile(dataDir, allFiles{idx}); end
    end
end

function ir = extractIR(S)
    ir = [];
    if isfield(S,'RiR') && ~isempty(S.RiR), ir = double(S.RiR(:));
    elseif isfield(S,'RIR') && ~isempty(S.RIR), ir = double(S.RIR(:));
    elseif isfield(S,'aufn') && ~isempty(S.aufn), ir = double(S.aufn(:));
    else
        fns = fieldnames(S);
        for f = 1:numel(fns)
            fname = fns{f};
            if startsWith(fname, '__'), continue; end
            v = S.(fname);
            if isnumeric(v) && numel(v) > 1000, ir = double(v(:)); return; end
        end
    end
end

function [ir_trunc, start_idx, end_idx, E_ratio, SNR_dB, dynamic_range_dB] = truncateIR(ir)
    ir_abs = abs(ir);
    max_amp = max(ir_abs);
    start_idx = find(ir_abs > max_amp * 0.01, 1, 'first'); % Start bei 1% vom Max (empfindlicher)
    if isempty(start_idx), start_idx = 1; end
    ir_trunc = ir(start_idx:end); % Nur Start abschneiden, Ende behalten für Decay
    end_idx = length(ir); E_ratio=0; SNR_dB=0; dynamic_range_dB=0;
end