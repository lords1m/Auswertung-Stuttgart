%% ============================================================
%  Darstellung_3D_Animation.m
%
%  Erstellt 3D-Animationen (Surface Plot) aus den Impulsantworten.
%  Z-Achse = Pegel in dB.
%
%  Output: .mp4 Dateien im Ordner 'Videos_3D'
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
outputVideoDir = 'Videos_3D'; 
fs = 500e3; % 500 kHz

% Video-Parameter
videoFPS = 20;          
timeStep_ms = 0.01;      % 0,1 Sekunden Schritte
windowSize_ms = 0.1;    
maxDuration_s = 0.1;    
cLim = [-60 0];         % Farbskala und Z-Achsen-Limit

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
    
    % 1. Daten laden
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
                    [ir_trunc, ~, ~, ~, ~, ~] = truncateIR(rawIR);
                    irs{r,c} = ir_trunc;
                catch, irs{r,c} = []; end
            else, irs{r,c} = []; end
        end
    end
    
    % 2. Zeitvektor
    dt = timeStep_ms / 1000;
    winLen = windowSize_ms / 1000;
    
    currentMaxLen = 0;
    for r=1:rows, for c=1:cols, currentMaxLen = max(currentMaxLen, length(irs{r,c})); end, end
    duration = min(maxDuration_s, currentMaxLen / fs);
    timePoints = 0:dt:duration;
    
    if isempty(timePoints), continue; end
    
    % 3. Video & 3D-Plot Initialisierung
    videoName = fullfile(outputVideoDir, ['3D_Animation_' variante '.mp4']);
    vObj = VideoWriter(videoName, 'MPEG-4');
    vObj.FrameRate = videoFPS;
    open(vObj);
    
    fig = figure('Visible', 'off', 'Position', [100, 100, 800, 600], 'Color', 'w');
    
    % Gitter für Surface Plot erstellen
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Initialer Surface Plot (flach bei -60dB)
    initialZ = ones(rows, cols) * cLim(1);
    hSurf = surf(X, Y, initialZ);
    
    % 3D-Optik Einstellungen
    shading interp;          % Glatte Farben (statt 'faceted' oder 'flat')
    colormap(jet);
    caxis(cLim);
    colorbar;
    
    % Achsen fixieren
    zlim([cLim(1), 10]);     % Z-Achse von -60 bis +10 dB
    xlim([1 cols]);
    ylim([1 rows]);
    
    % Ansicht einstellen (Azimuth, Elevation)
    view(-35, 45); 
    
    % Beschriftung
    xlabel('Spalten');
    ylabel('Zeilen');
    zlabel('Pegel [dBFS]');
    grid on;
    
    % Titel-Objekt
    hTitle = title('', 'FontSize', 14, 'FontWeight', 'bold');
    
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
                    val = cLim(1); % Stille = Minimum
                else
                    curr_idx_end = min(length(ir), idx_end);
                    segment = ir(idx_start:curr_idx_end);
                    rms_val = sqrt(mean(segment.^2));
                    val = 20 * log10((rms_val + eps) / MaxAmp_global);
                    if val < cLim(1), val = cLim(1); end
                end
                gridData(r,c) = val;
            end
        end
        
        % 3D-Daten aktualisieren
        set(hSurf, 'ZData', gridData, 'CData', gridData);
        set(hTitle, 'String', sprintf('%s\nZeit: %.3f s', strrep(variante,'_',' '), t_start));
        
        frame = getframe(fig);
        writeVideo(vObj, frame);
        
        % Fortschritt
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
    start_idx = find(ir_abs > max_amp * 0.01, 1, 'first'); 
    if isempty(start_idx), start_idx = 1; end
    ir_trunc = ir(start_idx:end); 
    end_idx = length(ir); E_ratio=0; SNR_dB=0; dynamic_range_dB=0;
end