%% ============================================================
%  Terzpegel-Auswertung (1/3 Oktav) aus RIR
%  Ergebnis: dBFS-Werte pro Terzband
%  Erweiterungen:
%   - automatische Referenzsuche (global über ausgewählte Varianten, _alt ausgeschlossen)
%   - selektive Auswahl von Varianten und Positionen
%   - Ausgabe: alle Plots in Plots/, alle Excels in Excel/
%   - Plot-Modus: 'absolute' oder 'difference' (Differenz gegen positionsmittleren Pegel)
% ============================================================

clear;
clc;

% Stelle sicher, dass wir im richtigen Verzeichnis sind
% (Falls Skript von anderem Ort aufgerufen wird)
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir)
    cd(scriptDir);
    fprintf('Arbeitsverzeichnis: %s\n', pwd);
else
    fprintf('Arbeitsverzeichnis: %s\n', pwd);
end

%% ---------------- Einstellungen (anpassen) ----------------
% Wenn empty -> automatisch alle Ordner 'Variante*_data' verwenden (ausser '*_alt')
selectedVariants = {};     % e.g. {'Variante_1_data','Variante_3_data'} oder {} für Auto
excludePattern = '*_alt';  % Ausschlussmuster

% Positionen (Messpunkte) auswählen
selectedPositions = 1:14;   % z.B. [1 3 5] oder 1:14

% Ausgabeordner
outputPlotDir = 'Plots';
outputExcelDir = 'Excel';

% Plot-Modus: 'absolute' (Pegel in dBFS relativ zur globalen Referenz)
%            'difference' (je Messung: Pegel minus Mean-Pegel über Positionen für diese Variante)
plotMode = 'absolute';   % oder 'difference'

% Sampling
fs = 500e3;   % 500 kHz Abtastrate

% Normgerechte 1/3-Oktav-Mitttenfrequenzen (IEC 61260)
f_terz = double([ ...
    560 625 710 800 900 1000 1120 1250 1400 1600 1800 2000 ...
    2240 2500 2800 3150 3550 4000 4500 5000 5600 6300 7100 8000 ...
    9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
    28000 31500 35500 40000 45000 50000 56000 63000 71000 80000 ...
    90000 100000 112000 126000 ]);

%% ---------------- Varianten-Entdeckung ----------------
% Alle passende Variante-Ordner finden
allDirs = dir();
allVariantFolders = allDirs([allDirs.isdir]);
% Filter nach Namen, die mit 'Variante' beginnen
variantMask = startsWith({allVariantFolders.name}, 'Variante');
allVariantFolders = allVariantFolders(variantMask);

% Filtern nach '*_data' falls vorhanden
dataMask = contains({allVariantFolders.name}, '_data');
dataFolders = allVariantFolders(dataMask);

if isempty(dataFolders)
    % Fallback: alle 'Variante_*' verwenden
    dataFolders = allVariantFolders;
end

variantNames = {dataFolders.name};

% Entferne Ordner, die dem Ausschlussmuster entsprechen ('_alt')
if ~isempty(excludePattern)
    toKeep = ~contains(variantNames, '_alt');  % einfache Ausschlussmechanik
    variantNames = variantNames(toKeep);
end

% Wenn user spezifische Varianten gesetzt hat -> filter auf diese Liste
if ~isempty(selectedVariants)
    variantNames = intersect(variantNames, selectedVariants, 'stable');
end

if isempty(variantNames)
    error('Keine Varianten zum Verarbeiten gefunden (nach Filter).');
end

fprintf('Verarbeite Varianten:\n');
for i=1:numel(variantNames), fprintf('  - %s\n', variantNames{i}); end

%% Debug: Überprüfe, welche Dateien in den Ordnern vorhanden sind
fprintf('\n--- DEBUG: Dateien in den Ordnern ---\n');
fprintf('Aktuelles Verzeichnis: %s\n', pwd);

% Teste was() um zu sehen, was existiert
all_items = dir;
fprintf('Alle Einträge im aktuellen Verzeichnis:\n');
for i = 1:numel(all_items)
    if all_items(i).isdir && contains(all_items(i).name, 'Variante')
        fprintf('  [DIR] %s\n', all_items(i).name);
        % Versuche, Unterordner zu lesen
        try
            subdir = dir(all_items(i).name);
            mat_count = sum(endsWith({subdir.name}, '.mat'));
            fprintf('    → %d .mat-Dateien\n', mat_count);
        catch
            fprintf('    → Fehler beim Lesen\n');
        end
    end
end
fprintf('---\n\n');

%% ---------------- Globale Referenz (FS) über alle ausgewählten Dateien finden ----------------
% Flexibles Laden: sucht alle .mat-Dateien mit 'Pos' + Positionsnummer

FS_global = 0;
found_any = false;
checkedFiles = {};

fprintf('\nSuche globale Referenz in %d Varianten...\n', numel(variantNames));

for vi = 1:numel(variantNames)
    variantName = variantNames{vi};
    fprintf('  Prüfe Variante: %s\n', variantName);

    % Verwende Java-API zum Auflisten von Dateien — zuverlässig in MATLAB auf macOS
    matFiles = {};
    try
        jDir = java.io.File(variantName);
        if jDir.exists() && jDir.isDirectory()
            jFiles = jDir.listFiles();
            for jf = 1:length(jFiles)
                jn = char(jFiles(jf).getName());
                if endsWith(jn, '.mat')
                    matFiles{end+1} = jn; %#ok<SAGROW>
                end
            end
        end
    catch
        matFiles = {};
    end

    if isempty(matFiles)
        fprintf('    Keine .mat-Dateien gefunden\n');
        continue;
    end
    fprintf('    Gefundene .mat-Dateien: %d\n', numel(matFiles));
    
    for pos = selectedPositions
        % Suche nach Datei mit 'Pos' und Positionsnummer
        posStr = num2str(pos);
        filename = '';
        for mf = 1:numel(matFiles)
            fname = matFiles{mf};
            if contains(fname, 'Pos') && contains(fname, posStr)
                filename = fullfile(variantName, fname);
                break;
            end
        end
        if isempty(filename)
            if pos == 1, fprintf('    Pos %d: nicht gefunden\n', pos); end
            continue;
        end
        
        checkedFiles{end+1} = filename; %#ok<SAGROW>
        try
            S = load(filename);
        catch ME
            if pos == 1, fprintf('    Pos %d: Fehler beim Laden - %s\n', pos, ME.message); end
            continue;
        end

        ir = [];
        if isfield(S,'RiR') && ~isempty(S.RiR)
            ir = double(S.RiR(:));
        elseif isfield(S,'RIR') && ~isempty(S.RIR)
            ir = double(S.RIR(:));
        else
            % Fallback: erstes numerisches Feld mit vielen Elementen
            fns = fieldnames(S);
            for f = 1:numel(fns)
                fname = fns{f};
                if startswith(fname, '__'), continue; end
                v = S.(fname);
                if isnumeric(v) && numel(v) > 1000
                    ir = double(v(:));
                    break;
                end
            end
        end

        if isempty(ir) || numel(ir) < 2
            if pos == 1, fprintf('    Pos %d: IR nicht geladen\n', pos); end
            continue;
        end

        N = numel(ir);
        IR_fft = fft(ir);
        H_mag = abs(IR_fft(1:floor(N/2)+1));
        FS_global = max(FS_global, max(H_mag));
        found_any = true;
        if pos == 1, fprintf('    Pos %d: FS_global = %g\n', pos, FS_global); end
    end
end

if ~found_any || FS_global == 0
    if isempty(checkedFiles)
        error('Keine geeigneten .mat-Dateien im erwarteten Format gefunden. Geprüfte Varianten: %s', strjoin(variantNames,', '));
    else
        error('Keine gültigen Impulsantworten gefunden; FS_Global = 0. Beispiel geprüft: %s', checkedFiles{1});
    end
end
fprintf('✓ Globale Referenz (FS) ermittelt: %g (Spektral-Maximum) aus %d Dateien\n', FS_global, numel(checkedFiles));

%% ---------------- Ensure output dirs ----------------
if ~exist(outputPlotDir,'dir'), mkdir(outputPlotDir); end
if ~exist(outputExcelDir,'dir'), mkdir(outputExcelDir); end

%% ---------------- Verarbeitung pro Variante ----------------
nFreq = length(f_terz);

for vi = 1:numel(variantNames)
    variantName = variantNames{vi};
    fprintf('\n========================================\n');
    fprintf('Verarbeite %s...\n', variantName);
    fprintf('========================================\n');

    % Speichergrößen: Positionen dynamisch
    posList = selectedPositions;
    nPos = numel(posList);

    L_dBFS = NaN(nPos, nFreq);
    H_all = cell(nPos,1);
    freq_all = cell(nPos,1);

    for pi = 1:nPos
        pos = posList(pi);

        % Verwende Java-API zum Auflisten von Dateien für diese Variante
        matFiles = {};
        try
            jDir = java.io.File(variantName);
            if jDir.exists() && jDir.isDirectory()
                jFiles = jDir.listFiles();
                for jf = 1:length(jFiles)
                    jn = char(jFiles(jf).getName());
                    if endsWith(jn, '.mat')
                        matFiles{end+1} = jn; %#ok<SAGROW>
                    end
                end
            end
        catch
            matFiles = {};
        end
        
        filename = '';
        posStr = num2str(pos);
        for mf = 1:numel(matFiles)
            fname = matFiles{mf};
            if contains(fname, 'Pos') && contains(fname, posStr)
                filename = fullfile(variantName, fname);
                break;
            end
        end
        if isempty(filename)
            warning('Keine Datei für Pos %d in %s gefunden. Position übersprungen.', pos, variantName);
            continue;
        end

        S = load(filename);
        ir = [];
        if isfield(S,'RiR') && ~isempty(S.RiR)
            ir = double(S.RiR(:));
        elseif isfield(S,'RIR') && ~isempty(S.RIR)
            ir = double(S.RIR(:));
        else
            % Fallback: erstes numerisches Feld mit vielen Elementen
            fns = fieldnames(S);
            for f = 1:numel(fns)
                fname = fns{f};
                if startswith(fname, '__'), continue; end
                v = S.(fname);
                if isnumeric(v) && numel(v) > 1000
                    ir = double(v(:));
                    break;
                end
            end
        end
        if isempty(ir) || numel(ir) < 2
            warning('Datei %s enthält keinen erkennbaren IR-Vektor. Übersprungen.', filename);
            continue;
        end
        N = length(ir);
        IR_fft = fft(ir);
        freq = (0:N-1) * fs / N;

        H_mag = abs(IR_fft(1:floor(N/2)+1));
        freq_pos = freq(1:floor(N/2)+1);

        % dBFS relativ zur globalen Referenz
        H_dBFS = 20*log10((H_mag + eps) / FS_global);

        % Terzband-Auswertung
        for k = 1:nFreq
            f1 = f_terz(k) / 2^(1/6);
            f2 = f_terz(k) * 2^(1/6);
            if f2 >= fs/2
                L_dBFS(pi,k) = NaN;
                continue;
            end
            idx = (freq >= f1) & (freq <= f2);
            if ~any(idx)
                L_dBFS(pi,k) = NaN;
                continue;
            end
            Power_band = abs(IR_fft(idx)).^2;
            total_energy = sum(Power_band) / N;
            rms_band = sqrt(total_energy);
            L_dBFS(pi,k) = 20*log10((rms_band + eps) / FS_global );
        end

        H_all{pi} = H_dBFS;
        freq_all{pi} = freq_pos;
    end

    % Mittelwert pro Terzband (über die ausgewählten Positionen)
    % Energieadditiv mitteln über Positionen (wenn Werte NaN -> omit)
    valid_mask = ~isnan(L_dBFS);
    if any(valid_mask(:))
        L_mean_dBFS = 10*log10( sum( 10.^(L_dBFS./10) .* valid_mask, 1, 'double' ) ./ max(1,sum(valid_mask,1)) );
    else
        L_mean_dBFS = NaN(1,nFreq);
    end

    % Optional: Differenzmodus -> Pegel relativ zum Positionsmittel (pro Variante)
    if strcmpi(plotMode,'difference')
        % Subtrahiere Mittelwert pro Terzband (L_mean_dBFS) von jeder Position
        L_plot = bsxfun(@minus, L_dBFS, L_mean_dBFS);
    else
        L_plot = L_dBFS;
    end

    %% ---------------- Excel-Export ----------------
    % Tabellen für diese Variante speichern
    rowNames = arrayfun(@(x) sprintf('Pos_%02d', x), posList, 'UniformOutput', false);
    colNames = arrayfun(@(f) sprintf('F%.0f', f), f_terz, 'UniformOutput', false);

    T_LdBFS = array2table(L_plot, 'RowNames', rowNames, 'VariableNames', colNames);

    excelFile = fullfile(outputExcelDir, sprintf('%s_Terzpegel_dBFS.xlsx', variantName));
    try
        writetable(T_LdBFS, excelFile, 'WriteRowNames', true);
        fprintf('Excel geschrieben: %s\n', excelFile);
    catch ME
        warning('Konnte Excel nicht schreiben: %s\n', ME.message);
    end

    %% ---------------- Plots ----------------
    variantPlotDir = fullfile(outputPlotDir, variantName);
    if ~exist(variantPlotDir,'dir'), mkdir(variantPlotDir); end

    % Y-Achsen-Grenzen (global für diese Variante)
    y_min_terz = min(L_plot(:), [], 'omitnan');
    y_max_terz = max(L_plot(:), [], 'omitnan');
    if isempty(y_min_terz) || isempty(y_max_terz) || isnan(y_min_terz) || isnan(y_max_terz)
        y_range_terz = [-100 0];
    else
        y_range_terz = [floor(y_min_terz/10)*10, ceil(y_max_terz/10)*10];
    end

    for pi = 1:nPos
        pos = posList(pi);
        fig = figure('Visible','off','Position',[100,100,1000,500]);
        stairs(f_terz, L_plot(pi,:), 'LineWidth',2, 'Color',[0 0.4470 0.7410]);
        grid on; set(gca,'XScale','log');
        xlabel('Frequenz [Hz]'); ylabel('Pegel [dB]');
        title(sprintf('%s - Position %02d (%s)', variantName, pos, plotMode));
        xlim([min(f_terz) max(f_terz)]);
        ylim(y_range_terz);
        xticks([500 1000 2000 5000 10000 20000 50000 100000]);
        xticklabels({'500','1k','2k','5k','10k','20k','50k','100k'});

        filename = fullfile(variantPlotDir, sprintf('Terzpegel_%s_Pos_%02d.png', variantName, pos));
        saveas(fig, filename);
        saveas(fig, strrep(filename,'.png','.fig'));
        close(fig);
    end

    % Übertragungsfunktionen für jede Position
    % Grenzen bestimmen
    y_min_H = inf; y_max_H = -inf;
    for pi = 1:nPos
        if ~isempty(H_all{pi})
            y_min_H = min(y_min_H, min(H_all{pi}));
            y_max_H = max(y_max_H, max(H_all{pi}));
        end
    end
    if ~isfinite(y_min_H), y_min_H = -120; end
    if ~isfinite(y_max_H), y_max_H = 0; end
    y_range_H = [floor(y_min_H/10)*10, ceil(y_max_H/10)*10];

    for pi = 1:nPos
        pos = posList(pi);
        if isempty(H_all{pi}), continue; end
        fig = figure('Visible','off','Position',[100,100,800,500]);
        semilogx(freq_all{pi}, H_all{pi}, '-', 'LineWidth',1.5, 'Color',[0.85 0.33 0.10]);
        grid on;
        xlabel('Frequenz [Hz]'); ylabel('Pegel [dBFS]');
        title(sprintf('Uebertragungsfunktion - %s - Pos %02d', variantName, pos));
        xlim([min(f_terz) max(f_terz)]);
        ylim(y_range_H);
        xticks([500 1000 2000 5000 10000 20000 50000 100000]);
        xticklabels({'500','1k','2k','5k','10k','20k','50k','100k'});

        filename = fullfile(variantPlotDir, sprintf('Uebertragungsfunktion_%s_Pos_%02d.png', variantName, pos));
        saveas(fig, filename);
        saveas(fig, strrep(filename,'.png','.fig'));
        close(fig);
    end

    fprintf('%s abgeschlossen: Excel + Plots gespeichert.\n', variantName);
end

disp('========================================');
disp('Alle ausgewählten Varianten erfolgreich verarbeitet!');
disp('========================================');


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:8540ce9e]
%   data: {"dataType":"text","outputData":{"text":"Terzpegel-Auswertung abgeschlossen.\n","truncated":false}}
%---
