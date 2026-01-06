%% VISUALISIERUNG: HIGH-SPEED PROPAGATION (0.0 - 0.1s)
% Fokus: Darstellung der Wellenfront in den ersten 100ms.
% Features: Interaktive Wahl von Variante & Terzband.

clear; clc; close all;

%% 1. KONFIGURATION
inputFolder = 'data';
outputFolder = 'Videos_HighSpeed_01s';
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

% Audio-Parameter
targetDuration = 0.1; % Wir laden EXAKT nur die ersten 0.1 Sekunden
fs_assumed = 48000;   % Fallback, falls fs nicht in Datei steht

% Varianten (Namen anpassen!)
variantNames = {'Variante_1', 'Variante_2', 'Variante_3', 'Variante_4'};

% Terzbänder (Mittenfrequenzen)
terzFreqs = [4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000, ...
             25000, 31500, 40000, 50000, 63000];

%% 2. INTERAKTIVE AUSWAHL
% Variante wählen
[vIdx, ok1] = listdlg('PromptString', 'Variante wählen:', ...
                      'SelectionMode', 'single', 'ListString', variantNames, ...
                      'ListSize', [200, 150]);
if ~ok1, return; end
selVar = variantNames{vIdx};

% Frequenz wählen
freqStr = string(terzFreqs) + " Hz";
[fIdx, ok2] = listdlg('PromptString', 'Terzband wählen:', ...
                      'SelectionMode', 'single', 'ListString', freqStr, ...
                      'ListSize', [200, 300]);
if ~ok2, return; end
selFreq = terzFreqs(fIdx);

%% 3. GEOMETRIE (Wand-Layout)
% Q1 = (0,0) unten links
cols_x = [0.5, 1.5, 2.5, 3.5]; 
rows_y = [3.0, 2.0, 1.0, 0.3]; 

% Mapping: [PosID, X, Y]
micMap = [
    1, cols_x(1), rows_y(1);   2, cols_x(2), rows_y(1);   3, cols_x(3), rows_y(1);   4, cols_x(4), rows_y(1);
    5, cols_x(1), rows_y(2);   6, cols_x(2), rows_y(2);   7, cols_x(3), rows_y(2);   8, cols_x(4), rows_y(2);
    9, cols_x(1), rows_y(3);  10, cols_x(2), rows_y(3);  11, cols_x(3), rows_y(3);  12, cols_x(4), rows_y(3);
   13, cols_x(2), rows_y(4);  14, cols_x(3), rows_y(4);  15, cols_x(4), rows_y(4);
];
numMics = size(micMap, 1);

%% 4. DATEN LADEN (Nur 0-100ms)
fprintf('Lade erste 0.1s für %s...\n', selVar);

% Wir checken erst die Samplerate der ersten Datei
filePattern = [selVar, '_Pos_%d.mat'];
testFile = fullfile(inputFolder, sprintf(filePattern, 1));

actualFs = fs_assumed;
if exist(testFile, 'file')
    info = load(testFile);
    if isfield(info, 'fs'), actualFs = info.fs; end
end

numSamples = round(targetDuration * actualFs);
rawData = zeros(numSamples, numMics);

for i = 1:numMics
    posID = micMap(i,1);
    fName = fullfile(inputFolder, sprintf(filePattern, posID));
    
    if exist(fName, 'file')
        tmp = load(fName);
        fn = fieldnames(tmp);
        sig = double(tmp.(fn{1}));
        
        % Nur die ersten numSamples nehmen
        L = min(length(sig), numSamples);
        if L > 0
            rawData(1:L, i) = sig(1:L);
        end
    else
        % Simulation (Backup)
        dist = sqrt(micMap(i,2)^2 + micMap(i,3)^2);
        delay = round((dist/343)*actualFs);
        if delay < numSamples
             rawData(delay, i) = 1; % Dirac
        end
    end
end

%% 5. FILTERUNG & ENVELOPE
fprintf('Filtere auf %d Hz...\n', selFreq);

% 5.1 Nyquist Check
if selFreq * 1.15 > actualFs/2
    error('Frequenz %d Hz ist zu hoch für Samplerate %d Hz (Nyquist)!', selFreq, actualFs);
end

% 5.2 Terz-Filter (Butterworth Bandpass)
f_lo = selFreq * 2^(-1/6);
f_hi = selFreq * 2^(1/6);
[b, a] = butter(4, [f_lo f_hi]/(actualFs/2), 'bandpass');

% Wichtig: filtfilt verwenden für 0 Phasenverschiebung (Zeitkorrektheit!)
filtData = filtfilt(b, a, rawData);

% 5.3 Hüllkurve (Envelope)
% Wir glätten das Signal extrem fein, um die Energiefront zu sehen.
% Fenster: ca. 0.2 ms (damit Impulse scharf bleiben)
winS = round(0.0002 * actualFs); 
if winS < 1, winS = 1; end

% Envelope: Wurzel aus gleitendem Mittelwert des Quadrats
envData = sqrt(movmean(filtData.^2, winS, 1));

% Normierung auf Max dieses Zeitfensters (für beste Sichtbarkeit)
envData = envData / (max(envData(:)) + 1e-9);

%% 6. VIDEO RENDERING (Zeitlupe)
fprintf('Rendere Video...\n');

% Grid Interpolation
pad = 0.5;
xLin = linspace(-pad, max