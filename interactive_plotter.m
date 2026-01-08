function interactive_plotter()
    % INTERACTIVE_PLOTTER Startet eine GUI zum Vergleichen von Messdaten
    % Liest Dateien aus dem Ordner 'processed' (Format: Proc_*.mat)
    
    % --- Konfiguration ---
    addpath('functions'); % Helper-Funktionen verfügbar machen
    procDir = 'processed';
    dataDir = 'data';
    if ~exist(procDir, 'dir')
        errordlg('Ordner "processed" nicht gefunden. Bitte erst step1 ausführen.', 'Fehler');
        return;
    end
    
    fileList = {}; % Wird durch updateFileList gefüllt

    isPlaying = false; % Status für Animation
    % Cache für Heatmap-Daten (um Neuladen beim Slider-Bewegen zu vermeiden)
    heatmapCache = struct('variante1', '', 'data1', [], 'variante2', '', 'data2', []); % data ist jetzt Map

    % Versuche, die globale Referenz (FS_global) aus einer verarbeiteten Datei zu laden,
    % damit die Raw-Daten im Spektrum identisch skaliert sind wie die Processed-Daten.
    FS_global_ref = 1.0;
    procFiles = dir(fullfile(procDir, 'Proc_*.mat'));
    if ~isempty(procFiles)
        try
            % Lade nur die Metadaten der ersten gefundenen Datei
            tmpLoad = load(fullfile(procDir, procFiles(1).name), 'Result');
            if isfield(tmpLoad.Result.meta, 'FS_global_used')
                FS_global_ref = tmpLoad.Result.meta.FS_global_used;
            end
        catch
        end
    end

    % --- GUI Layout erstellen ---
    f = figure('Name', 'Akustik Auswertung - Interaktiv', ...
               'NumberTitle', 'off', ...
               'Position', [100, 100, 1200, 700], ...
               'Color', 'w');

    % Panel für Steuerungen (Links)
    pnlControl = uipanel(f, 'Position', [0.01 0.01 0.25 0.98], 'Title', 'Einstellungen');

    % Modus Auswahl (nach oben verschoben)
    uicontrol(pnlControl, 'Style', 'text', 'Position', [10 660 200 20], 'String', 'Modus:', 'HorizontalAlignment', 'left');
    hMode = uicontrol(pnlControl, 'Style', 'popupmenu', 'Position', [10 635 200 25], ...
                      'String', {'Einzelansicht', 'Vergleich (Differenz)'}, ...
                      'Callback', @updateVisibility);

    % --- Messung 1 ---
    uicontrol(pnlControl, 'Style', 'text', 'Position', [10 605 200 20], 'String', 'Messung 1 (Referenz):', 'HorizontalAlignment', 'left');
    hSource1 = uicontrol(pnlControl, 'Style', 'popupmenu', 'Position', [10 580 200 25], ...
                      'String', {'Processed (Proc_*.mat)', 'Raw (Data/*.mat)'}, ...
                      'Callback', @updateFile1List);
    
    hFile1 = uicontrol(pnlControl, 'Style', 'listbox', 'Position', [10 430 200 140], ...
                       'String', fileList, 'Callback', @updatePlot);

    % --- Messung 2 ---
    lblFile2 = uicontrol(pnlControl, 'Style', 'text', 'Position', [10 400 200 20], 'String', 'Messung 2 (Vergleich):', 'HorizontalAlignment', 'left');
    
    hSource2 = uicontrol(pnlControl, 'Style', 'popupmenu', 'Position', [10 375 200 25], ...
                      'String', {'Processed (Proc_*.mat)', 'Raw (Data/*.mat)'}, ...
                      'Callback', @updateFile2List, 'Enable', 'off');
                      
    hFile2 = uicontrol(pnlControl, 'Style', 'listbox', 'Position', [10 225 200 140], ...
                       'String', fileList, 'Callback', @updatePlot, 'Enable', 'off');

    % Plot Typ
    uicontrol(pnlControl, 'Style', 'text', 'Position', [10 200 200 20], 'String', 'Darstellung:', 'HorizontalAlignment', 'left');
    hType = uicontrol(pnlControl, 'Style', 'popupmenu', 'Position', [10 180 200 25], ...
                      'String', {'Frequenzspektrum (Terz)', 'Impulsantwort (Zeit)', 'Energie über Zeit (ETC)', 'Energy Decay Curve (EDC)', 'Pegel über Entfernung', '3D Scatter (Raum)', 'Heatmap (Raumzeit)', 'Nachhallzeit (T30) über Frequenz'}, ...
                      'Callback', @updatePlot);

    % Frequenzfilter Checkbox
    hFilterFreq = uicontrol(pnlControl, 'Style', 'checkbox', 'Position', [10 140 200 20], ...
                            'String', 'Nur 4 kHz - 60 kHz', 'Value', 1, ...
                            'Callback', @updatePlot);

    % Feste Skalierung Checkbox
    hFixedScale = uicontrol(pnlControl, 'Style', 'checkbox', 'Position', [10 110 200 20], ...
                            'String', 'Feste Y-Achse (Global)', 'Value', 1, ...
                            'Callback', @updatePlot);

    % Feste X-Achse Checkbox
    hFixedX = uicontrol(pnlControl, 'Style', 'checkbox', 'Position', [10 80 200 20], ...
                            'String', 'Feste X-Achse (0-100ms)', 'Value', 0, ...
                            'Callback', @updatePlot);

    % Energie-Modus Checkbox (für Pegel über Entfernung)
    hEnergyMode = uicontrol(pnlControl, 'Style', 'checkbox', 'Position', [10 50 200 20], ...
                            'String', 'Energie (Linear) statt dB', 'Value', 0, ...
                            'Callback', @updatePlot);

    % --- Heatmap Controls (versteckt bis benötigt) ---
    hSliderTime = uicontrol(pnlControl, 'Style', 'slider', 'Position', [10 125 200 20], ...
                            'Min', -5, 'Max', 100, 'Value', 0, 'Visible', 'off', ...
                            'SliderStep', [0.001 0.05], ...
                            'Callback', @updateHeatmapFrame);
    lblTime = uicontrol(pnlControl, 'Style', 'text', 'Position', [10 145 200 20], ...
                        'String', 'Zeit: 0.0 ms', 'Visible', 'off', 'HorizontalAlignment', 'left');
    
    % Schwellenwert-Eingabe für Heatmap
    lblThreshold = uicontrol(pnlControl, 'Style', 'text', 'Position', [10 165 50 20], ...
                             'String', 'Min dB:', 'Visible', 'off', 'HorizontalAlignment', 'left');
    hEditThreshold = uicontrol(pnlControl, 'Style', 'edit', 'Position', [65 168 50 20], ...
                               'String', '-60', 'Visible', 'off', ...
                               'Callback', @updatePlot);

    hBtnPlay = uicontrol(pnlControl, 'Style', 'pushbutton', 'Position', [10 90 200 30], ...
                         'String', 'Play Animation', 'Visible', 'off', ...
                         'Callback', @playAnimation);

    % Speichern Button
    % Positionen angepasst, um Platz für Slider zu machen
    uicontrol(pnlControl, 'Style', 'pushbutton', 'Position', [10 50 200 30], ...
              'String', 'Plot speichern', 'Callback', @savePlot);

    % Batch Export Button
    uicontrol(pnlControl, 'Style', 'pushbutton', 'Position', [10 10 200 30], ...
              'String', 'Batch Export (Alle)', 'Callback', @batchExport);

    % Initialer Aufruf
    updateFile1List();
    updateFile2List();

    % --- Callbacks & Logik ---

    function updateFile1List(~, ~)
        updateListGeneric(hSource1, hFile1);
        updatePlot();
    end

    function updateFile2List(~, ~)
        updateListGeneric(hSource2, hFile2);
        updatePlot();
    end

    function updateListGeneric(hSrc, hList)
        if hSrc.Value == 1
            searchDir = procDir;
            pattern = 'Proc_*.mat';
        else
            searchDir = dataDir;
            pattern = '*.mat';
        end
        
        newList = {};
        if ~exist(searchDir, 'dir')
             newList = {'Ordner nicht gefunden'};
        else
             files = dir(fullfile(searchDir, pattern));
             if isempty(files)
                 newList = {'Keine Dateien'};
             else
                 newList = {files.name};
             end
        end
        set(hList, 'String', newList, 'Value', 1);
    end

    function updateVisibility(~, ~)
        mode = hMode.Value;
        if mode == 1 % Einzel
            set(hFile2, 'Enable', 'off');
            set(hSource2, 'Enable', 'off');
            set(lblFile2, 'Enable', 'off');
        else % Vergleich
            set(hFile2, 'Enable', 'on');
            set(hSource2, 'Enable', 'on');
            set(lblFile2, 'Enable', 'on');
        end
        updatePlot();
    end

    function updatePlot(~, ~)
        % Listen holen
        list1 = get(hFile1, 'String');
        list2 = get(hFile2, 'String');
        
        if isempty(list1) || strcmp(list1{1}, 'Keine Dateien'), cla(findobj(f, 'Type', 'axes')); return; end

        % 1. Daten laden
        idx1 = hFile1.Value;
        name1 = list1{idx1};
        R1 = loadData(name1, hSource1.Value);

        isCompare = (hMode.Value == 2);
        R2 = [];
        if isCompare
            if isempty(list2) || strcmp(list2{1}, 'Keine Dateien'), return; end
            idx2 = hFile2.Value;
            name2 = list2{idx2};
            R2 = loadData(name2, hSource2.Value);
        end

        % ... Rest der Funktion bleibt gleich bis zu den Limits ...
        
        plotType = hType.Value; % 1=Spec, 2=IR, 3=ETC
        useFreqFilter = hFilterFreq.Value;
        useFixedScale = hFixedScale.Value;
        useFixedX = hFixedX.Value;
        useEnergyMode = hEnergyMode.Value;

        % UI Sichtbarkeit steuern
        if plotType == 7 % Heatmap
            set(hSliderTime, 'Visible', 'on');
            set(lblTime, 'Visible', 'on');
            set(hBtnPlay, 'Visible', 'on');
            set(lblThreshold, 'Visible', 'on');
            set(hEditThreshold, 'Visible', 'on');
            % Verstecke andere Controls die stören könnten
            set(hFixedScale, 'Visible', 'off');
            set(hFixedX, 'Visible', 'off');
            set(hFilterFreq, 'Visible', 'off'); % Filter stört hier
        else
            set(hSliderTime, 'Visible', 'off'); set(lblTime, 'Visible', 'off'); set(hBtnPlay, 'Visible', 'off');
            set(lblThreshold, 'Visible', 'off');
            set(hEditThreshold, 'Visible', 'off');
            set(hFixedScale, 'Visible', 'on');
            set(hFilterFreq, 'Visible', 'on');
        end

        % --- AXES MANAGEMENT ---
        % Alte Achsen löschen, um Konflikte (Subplot vs. Single) zu vermeiden
        delete(findobj(f, 'Type', 'axes'));

        if plotType == 1 && isCompare
            % Layout für Vergleich: Zwei Plots übereinander
            ax1 = axes(f, 'Position', [0.32 0.55 0.65 0.38]); 
            ax2 = axes(f, 'Position', [0.32 0.10 0.65 0.35]);
        else
            % Layout Einzel: Ein großer Plot
            ax = axes(f, 'Position', [0.32 0.1 0.65 0.85]);
            grid(ax, 'on'); hold(ax, 'on');
        end

        % --- PLOT LOGIK ---
        switch plotType
            case 1 % SPEKTRUM (Terz)
                f_vec = R1.freq.f_center;
                y1 = R1.freq.terz_dbfs;
                
                % Filter anwenden (4k - 60k)
                if useFreqFilter
                    mask = (f_vec >= 4000) & (f_vec <= 60000);
                else
                    mask = true(size(f_vec));
                end
                
                f_sub = f_vec(mask);
                y1_sub = y1(mask);
                x_idx = 1:length(f_sub); % Indizes für gleichmäßige Darstellung
                
                % Labels für X-Achse generieren
                x_labels = arrayfun(@(x) sprintf('%g', x), f_sub, 'UniformOutput', false);
                
                if isCompare
                    y2 = R2.freq.terz_dbfs;
                    y2_sub = y2(mask);
                    
                    % Oben: Beide Spektren
                    axes(ax1); 
                    name1_leg = sprintf('%s (L_{sum}=%.1f dB)', cleanName(name1), R1.freq.sum_level);
                    name2_leg = sprintf('%s (L_{sum}=%.1f dB)', cleanName(name2), R2.freq.sum_level);
                    
                    stairs(x_idx, y1_sub, 'b-', 'LineWidth', 1.5, 'DisplayName', name1_leg); hold on;
                    stairs(x_idx, y2_sub, 'r-', 'LineWidth', 1.5, 'DisplayName', name2_leg);
                    grid on; legend show; ylabel('Pegel [dBFS]'); title('Frequenzgang Vergleich');
                    set(gca, 'XTick', x_idx, 'XTickLabel', x_labels);
                    xtickangle(45);
                    xlim([0.5 length(f_sub)+0.5]);
                    if useFixedScale
                        ylim([-30 5]); % Fester Bereich angepasst auf -30dB
                    end
                    
                    % Unten: Differenz
                    axes(ax2);
                    diff_y = y2_sub - y1_sub;
                    bar(x_idx, diff_y, 'FaceColor', [0.5 0.5 0.5]);
                    grid on; ylabel('Differenz [dB]'); xlabel('Frequenz [Hz]');
                    title('Differenz (Messung 2 - Messung 1)');
                    set(gca, 'XTick', x_idx, 'XTickLabel', x_labels);
                    xtickangle(45);
                    xlim([0.5 length(f_sub)+0.5]);
                else
                    % Einzelplot
                    % ax ist bereits aktiv
                    stairs(x_idx, y1_sub, 'b-', 'LineWidth', 2);
                    grid on; xlabel('Frequenz [Hz]'); ylabel('Pegel [dBFS]');
                    title(sprintf('Spektrum: %s (L_{sum} = %.1f dB)', cleanName(name1), R1.freq.sum_level));
                    set(gca, 'XTick', x_idx, 'XTickLabel', x_labels);
                    xtickangle(45);
                    xlim([0.5 length(f_sub)+0.5]);
                    if useFixedScale
                        ylim([-30 5]);
                    end
                end

            case 2 % IMPULSANTWORT (Zeit)
                t1 = (0:length(R1.time.ir)-1) / R1.meta.fs * 1000; % ms
                
                if isCompare
                    t2 = (0:length(R2.time.ir)-1) / R2.meta.fs * 1000;
                    plot(t1, R1.time.ir, 'b', 'DisplayName', cleanName(name1));
                    plot(t2, R2.time.ir, 'r', 'DisplayName', cleanName(name2));
                    legend show;
                    title('Impulsantworten Vergleich');
                else
                    if all(R1.time.ir == 0)
                        text(0.5, 0.5, 'Keine Zeitdaten (Average File)', 'Units', 'normalized', 'HorizontalAlignment', 'center');
                    else
                        plot(t1, R1.time.ir, 'b');
                    end
                    title(['Impulsantwort: ' cleanName(name1)]);
                    % Info über Truncation anzeigen
                    infoStr = sprintf('Truncation: Idx %d bis %d\nEnergieanteil: %.1f%%', ...
                        R1.time.metrics.idx_start, R1.time.metrics.idx_end, R1.time.metrics.energy_share*100);
                    text(0.05, 0.9, infoStr, 'Units', 'normalized', 'BackgroundColor', 'w');
                end
                xlabel('Zeit [ms]'); ylabel('Amplitude');
                grid on;
                if useFixedScale
                    ylim([-FS_global_ref, FS_global_ref] * 1.1);
                end
                if useFixedX
                    xlim([0 100]);
                end

            case 3 % ENERGIE ÜBER ZEIT (ETC)
                % Berechnung der Hüllkurve (Logarithmisch)
                % Wir nutzen hier vereinfacht den Betrag in dB, geglättet wäre besser, 
                % aber für Rohdaten-Check reicht das.
                
                etc1 = 20*log10(abs(R1.time.ir) + eps);
                t1 = (0:length(etc1)-1) / R1.meta.fs * 1000;
                
                % Dynamikbereich begrenzen für Plot (z.B. max - 60dB)
                maxVal1 = max(etc1);
                
                if isCompare
                    etc2 = 20*log10(abs(R2.time.ir) + eps);
                    t2 = (0:length(etc2)-1) / R2.meta.fs * 1000;
                    maxVal2 = max(etc2);
                    globalMax = max(maxVal1, maxVal2);
                    
                    plot(t1, etc1, 'b', 'DisplayName', cleanName(name1));
                    plot(t2, etc2, 'r', 'DisplayName', cleanName(name2));
                    legend show;
                    ylim([globalMax-60, globalMax+5]);
                    title('Energie-Zeit-Kurve (ETC) Vergleich');
                else
                    plot(t1, etc1, 'b');
                    ylim([maxVal1-60, maxVal1+5]);
                    title(['Energie-Zeit-Kurve (ETC): ' cleanName(name1)]);
                end
                xlabel('Zeit [ms]'); ylabel('Pegel [dB]');
                grid on;
                if useFixedScale
                    ylim([-30 5]);
                end
                if useFixedX
                    xlim([0 100]);
                end

            case 4 % ENERGY DECAY CURVE (EDC)
                edc1 = calc_edc(R1.time.ir);
                t1 = (0:length(edc1)-1) / R1.meta.fs * 1000;
                
                if isCompare
                    edc2 = calc_edc(R2.time.ir);
                    t2 = (0:length(edc2)-1) / R2.meta.fs * 1000;
                    
                    plot(t1, edc1, 'b', 'DisplayName', cleanName(name1));
                    plot(t2, edc2, 'r', 'DisplayName', cleanName(name2));
                    legend show;
                    title('Energy Decay Curve (EDC) Vergleich');
                else
                    plot(t1, edc1, 'b');
                    title(['Energy Decay Curve (EDC): ' cleanName(name1)]);
                end
                xlabel('Zeit [ms]'); ylabel('Pegel [dB]');
                if useFixedScale
                    ylim([-30 5]); 
                end
                if useFixedX
                    xlim([0 100]);
                end
                grid on;

            case 5 % PEGEL ÜBER ENTFERNUNG
                % Daten sammeln für Variante 1 (basiert auf der Variante der ausgewählten Datei)
                [dist1, lev1, ~] = get_variant_levels(R1.meta.variante, hSource1.Value, useEnergyMode);
                
                % Plot 1
                scatter(dist1, lev1, 60, 'b', 'filled', 'DisplayName', strrep(R1.meta.variante, '_', ' ')); hold on;
                
                % Ideal Kurve berechnen (basierend auf Var 1)
                if ~isempty(dist1)
                    [dist_ideal, L_ideal] = calc_ideal_curve(dist1, lev1, useEnergyMode);
                    plot(dist_ideal, L_ideal, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal 1/r');
                end
                
                if isCompare
                    [dist2, lev2, ~] = get_variant_levels(R2.meta.variante, hSource2.Value, useEnergyMode);
                    scatter(dist2, lev2, 60, 'r', 'filled', 'DisplayName', strrep(R2.meta.variante, '_', ' '));
                end
                
                xlabel('Entfernung von Quelle [m]'); 
                if useEnergyMode
                    ylabel('Energie (linear)');
                    title('Energieabfall über Entfernung');
                else
                    ylabel('Summenpegel [dBFS]');
                    title('Pegelabfall über Entfernung');
                end
                grid on; legend show;
                
                if useFixedScale && ~useEnergyMode
                    ylim([-30 5]);
                end
                
                % Text-Info
                text(0.02, 0.05, 'Punkte zeigen Messpositionen 1-15', 'Units', 'normalized', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);

            case 6 % 3D SCATTER (RAUM)
                geo = get_geometry();
                
                % --- Daten 1 ---
                [dist1, lev1, pos1] = get_variant_levels(R1.meta.variante, hSource1.Value, useEnergyMode);
                x1 = []; y1 = [];
                for p = pos1
                    g = geo([geo.pos] == p);
                    if ~isempty(g), x1(end+1) = g.x; y1(end+1) = g.y; end
                end
                
                scatter3(x1, y1, lev1, 80, 'b', 'filled', 'DisplayName', strrep(R1.meta.variante, '_', ' ')); hold on;
                
                % Text Labels
                for i = 1:length(x1)
                    text(x1(i), y1(i), lev1(i), sprintf('  P%d', pos1(i)), 'FontSize', 8);
                end

                % --- Daten 2 (Vergleich) ---
                if isCompare
                    [dist2, lev2, pos2] = get_variant_levels(R2.meta.variante, hSource2.Value, useEnergyMode);
                    x2 = []; y2 = [];
                    for p = pos2
                        g = geo([geo.pos] == p);
                        if ~isempty(g), x2(end+1) = g.x; y2(end+1) = g.y; end
                    end
                    scatter3(x2, y2, lev2, 80, 'r', 'filled', 'DisplayName', strrep(R2.meta.variante, '_', ' '));
                    
                    % Verbindungslinien zwischen gleichen Positionen
                    common = intersect(pos1, pos2);
                    for p = common
                        i1 = find(pos1 == p);
                        i2 = find(pos2 == p);
                        plot3([x1(i1) x1(i1)], [y1(i1) y1(i1)], [lev1(i1) lev2(i2)], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
                    end
                end
                
                % Quelle visualisieren
                z_source = max(lev1); 
                if isCompare && ~isempty(lev2), z_source = max(z_source, max(lev2)); end
                scatter3(0, 0, z_source, 100, 'k', 'filled', 'DisplayName', 'Quelle');
                
                grid on; view(45, 30);
                xlabel('X [m]'); ylabel('Y [m]');
                if useEnergyMode
                    zlabel('Energie (linear)');
                    title('Räumliche Energieverteilung');
                else
                    zlabel('Pegel [dBFS]');
                    title('Räumliche Pegelverteilung');
                end
                legend show;
                
                if useFixedScale && ~useEnergyMode
                    zlim([-30 5]);
                end

            case 7 % HEATMAP (RAUMZEIT)
                t_ms = get(hSliderTime, 'Value');
                set(lblTime, 'String', sprintf('Zeit: %.1f ms', t_ms));
                
                % Daten laden (gecached)
                data1 = get_heatmap_data(R1.meta.variante, 1);
                grid1 = calc_heatmap_grid(data1, t_ms, R1.meta.fs);
                
                % Threshold aus GUI lesen
                min_db = str2double(get(hEditThreshold, 'String'));
                if isnan(min_db), min_db = -60; end
                cLim = [min_db 0]; 
                
                if isCompare
                    data2 = get_heatmap_data(R2.meta.variante, 2);
                    grid2 = calc_heatmap_grid(data2, t_ms, R2.meta.fs);
                    
                    % Plot 1
                    axes(ax1);
                    imagesc(grid1);
                    colormap(ax1, jet); caxis(ax1, cLim); colorbar;
                    title(sprintf('%s (%.1f ms)', strrep(R1.meta.variante,'_',' '), t_ms));
                    axis square; axis off;
                    add_heatmap_labels(grid1, min_db);
                    
                    % Plot 2
                    axes(ax2);
                    imagesc(grid2);
                    colormap(ax2, jet); caxis(ax2, cLim); colorbar;
                    title(sprintf('%s (%.1f ms)', strrep(R2.meta.variante,'_',' '), t_ms));
                    axis square; axis off;
                    add_heatmap_labels(grid2, min_db);
                else
                    % Einzelplot
                    imagesc(grid1);
                    colormap(jet); caxis(cLim); colorbar;
                    title(sprintf('Energieverteilung: %s @ %.1f ms', strrep(R1.meta.variante,'_',' '), t_ms));
                    axis square; axis off;
                    add_heatmap_labels(grid1, min_db);
                end

            case 8 % NACHHALLZEIT (T30) ÜBER FREQUENZ
                % Messung 1: Prüfen ob vorberechnet oder neu berechnen
                if isfield(R1.freq, 't30') && ~isempty(R1.freq.t30)
                    t30_1 = R1.freq.t30;
                    f_vec = R1.freq.t30_freqs;
                else
                    [t30_1, f_vec] = calc_rt60_spectrum(R1.time.ir, R1.meta.fs);
                end
                
                if isCompare
                    if isfield(R2.freq, 't30') && ~isempty(R2.freq.t30), t30_2 = R2.freq.t30;
                    else, [t30_2, ~] = calc_rt60_spectrum(R2.time.ir, R2.meta.fs); end
                    
                    plot(f_vec, t30_1, 'b-o', 'LineWidth', 1.5, 'DisplayName', cleanName(name1)); hold on;
                    plot(f_vec, t30_2, 'r-x', 'LineWidth', 1.5, 'DisplayName', cleanName(name2));
                    legend show;
                    title('Nachhallzeit T30 Vergleich');
                else
                    plot(f_vec, t30_1, 'b-o', 'LineWidth', 1.5);
                    title(['Nachhallzeit T30: ' cleanName(name1)]);
                end
                
                set(gca, 'XScale', 'log');
                xlabel('Frequenz [Hz]'); ylabel('Nachhallzeit T30 [s]');
                grid on; xlim([3500 65000]);
                if useFixedScale, ylim([0 2.5]); end
        end
    end

    function R = loadData(filename, sourceType)
        if sourceType == 1
            % Processed
            tmp = load(fullfile(procDir, filename));
            R = tmp.Result;
        else
            % Raw Data
            filepath = fullfile(dataDir, filename);
            [S, meta] = load_and_parse_file(filepath);
            ir = extract_ir(S);
            if isempty(ir), ir = zeros(100,1); end
            ir = ir - mean(ir); % DC-Offset entfernen für konsistente Darstellung
            
            R.time.ir = ir;
            R.time.metrics.idx_start = 1;
            R.time.metrics.idx_end = length(ir);
            R.time.metrics.energy_share = 1;
            R.meta = meta;
            if ~isfield(R.meta, 'fs'), R.meta.fs = 500e3; end
            
            % On-the-fly Berechnung des Spektrums (mit globaler Referenz für korrekte Skalierung)
            [L_terz, L_sum, f_center] = calc_terz_spectrum(R.time.ir, R.meta.fs, FS_global_ref);
            R.freq.f_center = f_center;
            R.freq.terz_dbfs = L_terz;
            R.freq.sum_level = L_sum;
        end
    end

    function [t60_vals, f_center] = calc_rt60_spectrum(ir, fs)
        % Frequenzen (4k - 60k)
        f_center = [4000 5000 6300 7100 8000 9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 28000 31500 35500 40000 45000 50000 56000];
        t60_vals = NaN(size(f_center));
        
        for k = 1:length(f_center)
            fc = f_center(k);
            % Bandgrenzen für Terzfilter
            fl = fc / 2^(1/6);
            fu = fc * 2^(1/6);
            
            if fu >= fs/2, continue; end
            
            try
                % Butterworth Bandpass (n=4 -> Ordnung 8)
                [b, a] = butter(4, [fl fu]/(fs/2), 'bandpass');
                filt_ir = filtfilt(b, a, ir);
                
                % EDC berechnen (lokale Funktion nutzen)
                edc_db = calc_edc(filt_ir);
                
                % T30 Bestimmung (-5 dB bis -35 dB)
                idx_start = find(edc_db <= -5, 1, 'first');
                
                if ~isempty(idx_start)
                    % Suche nach -35 dB nach dem Startpunkt
                    idx_end_rel = find(edc_db(idx_start:end) <= -35, 1, 'first');
                    
                    if ~isempty(idx_end_rel)
                        idx_end = idx_start + idx_end_rel - 1;
                        
                    % Lineare Regression auf dem Abschnitt
                    y_segment = edc_db(idx_start:idx_end);
                    x_segment = (0:length(y_segment)-1)' / fs; 
                    
                    if length(y_segment) > 5
                        p = polyfit(x_segment, y_segment, 1);
                        slope = p(1); % Steigung in dB/s
                        
                        if slope < 0, t60_vals(k) = -60 / slope; end
                    end
                    end
                end
            catch
            end
        end
    end
    
    function edc_db = calc_edc(ir)
        % Schroeder Integration (Backward Integration)
        % Energie berechnen
        E = cumsum(ir(end:-1:1).^2);
        E = E(end:-1:1);
        
        % Normalisieren auf 0 dB
        max_E = max(E);
        if max_E == 0
            edc_db = ones(size(ir)) * -100;
        else
            edc_db = 10*log10(E / max_E + eps);
        end
    end

    function [dist, levels, positions] = get_variant_levels(variante, sourceType, energyMode)
        geo = get_geometry();
        dist = []; levels = []; positions = [];
        
        % Dateiliste holen
        if sourceType == 1 % Processed
            files = dir(fullfile(procDir, sprintf('Proc_%s_Pos*.mat', variante)));
        else % Raw
            % Versuch, Dateien vorzufiltern (Wildcards helfen bei Performance)
            files = dir(fullfile(dataDir, sprintf('*%s*.mat', variante)));
        end
        
        for i = 1:length(files)
            try
                if sourceType == 1
                    D = load(fullfile(files(i).folder, files(i).name), 'Result');
                    meta = D.Result.meta;
                    if energyMode
                        val = D.Result.time.metrics.energy;
                    else
                        val = D.Result.freq.sum_level;
                    end
                else
                    [S, meta] = load_and_parse_file(fullfile(files(i).folder, files(i).name));
                    % Sicherstellen, dass es wirklich die richtige Variante ist
                    if ~strcmp(meta.variante, variante), continue; end
                    
                    ir = extract_ir(S);
                    if isempty(ir), continue; end
                    % Pegel berechnen (RMS über alles, DC-bereinigt)
                    if energyMode
                        val = sum((ir - mean(ir)).^2);
                    else
                        val = 20*log10(rms(ir - mean(ir)) / FS_global_ref + eps);
                    end
                end
                
                % Position parsen
                posNum = str2double(meta.position);
                if isnan(posNum), continue; end
                
                % Geometrie finden
                idx = find([geo.pos] == posNum);
                if ~isempty(idx)
                    d = geo(idx).distance;
                    if d > 0 % Quelle (d=0) ausschließen für Log-Plot
                        dist(end+1) = d;
                        levels(end+1) = val;
                        positions(end+1) = posNum;
                    end
                end
            catch
            end
        end
    end

    function [d_line, L_line] = calc_ideal_curve(dist, levels, energyMode)
        % Referenz: Mittelwert der kleinsten Distanz
        if isempty(dist), d_line=[]; L_line=[]; return; end
        
        min_d = min(dist);
        % Toleranzbereich für Referenzpunkte (z.B. alle bei 0.3m)
        ref_mask = abs(dist - min_d) < 0.05;
        L_ref = mean(levels(ref_mask));
        
        d_line = linspace(min_d, max(dist)*1.1, 100);
        
        if energyMode
            % E(r) ~ 1/r^2
            L_line = L_ref * (min_d ./ d_line).^2;
        else
            % L(r) = L_ref - 20*log10(r/r_ref) (Halbraum-Dämpfung)
            L_line = L_ref - 20*log10(d_line / min_d);
        end
    end

    function pos_info = get_geometry()
        % Hardcoded Geometrie (Grid 4x4, Quelle bei 0,0)
        pos_info = struct('pos', {}, 'x', {}, 'y', {}, 'distance', {});
        
        % Koordinaten (PosID, x, y)
        coords = [
            1, 0, 1.2; 2, 0.3, 1.2; 3, 0.6, 1.2; 4, 1.2, 1.2;
            5, 0, 0.6; 6, 0.3, 0.6; 7, 0.6, 0.6; 8, 1.2, 0.6;
            9, 0, 0.3; 10, 0.3, 0.3; 11, 0.6, 0.3; 12, 1.2, 0.3;
            13, 0.3, 0; 14, 0.6, 0; 15, 1.2, 0
        ];
        
        source_x = 0; source_y = 0;
        
        for i = 1:size(coords, 1)
            p = coords(i, 1);
            x = coords(i, 2);
            y = coords(i, 3);
            d = sqrt((x - source_x)^2 + (y - source_y)^2);
            
            pos_info(i).pos = p;
            pos_info(i).x = x;
            pos_info(i).y = y;
            pos_info(i).distance = d;
        end
    end

    function savePlot(~, ~)
        [filename, pathname] = uiputfile({'*.png', 'PNG Bild (*.png)'; '*.pdf', 'PDF Dokument (*.pdf)'; '*.fig', 'MATLAB Figur (*.fig)'}, 'Plot speichern als...');
        if isequal(filename, 0), return; end
        savePath = fullfile(pathname, filename);
        
        % Temporäre Figur für sauberen Export (ohne GUI-Elemente)
        f_temp = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1000 700]);
        new_ax = copyobj(findobj(f, 'Type', 'axes'), f_temp);
        
        % Colormap übertragen (wichtig für Heatmap)
        colormap(f_temp, colormap(f));
        
        % Positionen der Achsen anpassen (im GUI rechts, hier zentriert)
        for i = 1:length(new_ax)
            try
                pos = get(new_ax(i), 'Position');
                new_x = max(0.05, pos(1) - 0.22); % Verschiebe nach links statt Breite zu erzwingen
                set(new_ax(i), 'Position', [new_x, pos(2), pos(3), pos(4)]);
            catch
            end
        end
        
        try
            if endsWith(filename, '.fig')
                savefig(f_temp, savePath);
            elseif exist('exportgraphics', 'file')
                exportgraphics(f_temp, savePath, 'Resolution', 300);
            else
                saveas(f_temp, savePath);
            end
        catch
        end
        delete(f_temp);
    end

    function batchExport(~, ~)
        % 1. Ordner wählen
        outDir = uigetdir(pwd, 'Ordner für Batch-Export wählen');
        if isequal(outDir, 0), return; end
        
        % 2. Check
        list = get(hFile1, 'String');
        if ischar(list), list = {list}; end % Sicherstellen, dass list ein Cell-Array ist
        
        if isempty(list) || (length(list)==1 && strcmp(list{1}, 'Keine Dateien'))
            msgbox('Keine Dateien vorhanden.', 'Info'); return;
        end
        
        % Modus Warnung (Batch ist für Einzelansicht gedacht)
        if hMode.Value == 2
            answer = questdlg('Batch-Export wechselt zur Einzelansicht. Fortfahren?', 'Modus', 'Ja', 'Nein', 'Ja');
            if ~strcmp(answer, 'Ja'), return; end
            hMode.Value = 1;
            updateVisibility([], []);
        end
        
        % Export-Modus bestimmen: Variantenweise (für globale Plots) oder Dateiweise
        plotType = hType.Value;
        % Globale Plots: 5 (Pegel), 6 (3D), 7 (Heatmap)
        isGlobalPlot = ismember(plotType, [5, 6, 7]);
        
        if isGlobalPlot
            % Varianten identifizieren (nur einmal pro Variante exportieren)
            variants = {};
            indices = [];
            
            for i = 1:length(list)
                fname = list{i};
                % Bereinigen für Gruppierung
                cleanFname = regexprep(fname, '^Proc_', '');
                % Regex für Variantenname (alles vor _Pos oder _Quelle)
                tokens = regexp(cleanFname, '^(.*?)(?:[_,]Pos|[_,]Quelle|Pos|Quelle)', 'tokens', 'once', 'ignorecase');
                if ~isempty(tokens)
                    vName = tokens{1};
                else
                    [~, vName, ~] = fileparts(cleanFname);
                end
                
                if ~ismember(vName, variants)
                    variants{end+1} = vName;
                    indices(end+1) = i;
                end
            end
            
            exportList = variants;
            exportIndices = indices;
            modeLabel = 'Varianten';
        else
            exportList = list;
            exportIndices = 1:length(list);
            modeLabel = 'Dateien';
        end
        
        nExport = length(exportIndices);
        if nExport == 0
            msgbox('Nichts zu exportieren.', 'Info'); return;
        end
        
        hWait = waitbar(0, sprintf('Starte Export (%s)...', modeLabel));
        fprintf('Starte Batch-Export von %d Einträgen (%s)...\n', nExport, modeLabel);
        if isGlobalPlot
            fprintf('Hinweis: Globale Plots (Heatmap, 3D) fassen alle Positionen einer Variante zusammen.\n');
            fprintf('Gefundene Varianten: %d (aus %d Dateien)\n', length(variants), length(list));
            fprintf('Es wird daher nur EIN Bild pro Variante exportiert.\n');
        end
        
        try
            for k = 1:nExport
                if ~ishandle(hWait)
                    fprintf('Abbruch durch Benutzer.\n');
                    break; 
                end 
                
                idx = exportIndices(k);
                itemLabel = exportList{k};
                
                try
                    % Waitbar Update (sicher gegen Sonderzeichen)
                    cleanLabel = strrep(itemLabel, '_', '\_');
                    msg = sprintf('Exportiere %d/%d: %s', k, nExport, cleanLabel);
                    waitbar((k-1)/nExport, hWait, msg);
                catch
                    % Ignoriere Waitbar-Fehler
                end
                
                try
                    % Auswählen & Plotten
                    set(hFile1, 'Value', idx);
                    updatePlot([], []);
                    drawnow; 
                    pause(0.1); % Kurz warten für Rendering
                    
                    % Dateiname generieren
                    if isGlobalPlot
                        % Präfix basierend auf Plot-Typ
                        plotStrings = get(hType, 'String');
                        pName = plotStrings{plotType};
                        if contains(pName, 'Heatmap'), prefix = 'Heatmap';
                        elseif contains(pName, '3D'), prefix = '3D_Raum';
                        elseif contains(pName, 'Pegel'), prefix = 'Pegel_Entfernung';
                        else, prefix = 'Plot'; end
                        
                        % Dateiname bereinigen
                        safeLabel = regexprep(itemLabel, '[^a-zA-Z0-9_\-]', '');
                        baseName = fullfile(outDir, sprintf('%s_%s', prefix, safeLabel));
                    else
                        [~, fname, ~] = fileparts(list{idx});
                        baseName = fullfile(outDir, fname);
                    end
                    
                    fprintf('  [%d/%d] Speichere: %s ... ', k, nExport, baseName);
                    
                    % Export (Kopie der Achsen für sauberes Layout)
                    f_temp = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1000 700]);
                    
                    % Achsen kopieren
                    ax_src = findobj(f, 'Type', 'axes');
                    if isempty(ax_src)
                        fprintf('SKIPPED (Kein Plot)\n');
                        close(f_temp);
                        continue;
                    end
                    new_ax = copyobj(ax_src, f_temp);
                    
                    % Colormap übertragen
                    colormap(f_temp, colormap(f));
                    
                    for ax_i = 1:length(new_ax)
                        try
                            pos = get(new_ax(ax_i), 'Position');
                            new_x = max(0.05, pos(1) - 0.22); % Verschiebe nach links, erhalte Aspect Ratio
                            set(new_ax(ax_i), 'Position', [new_x, pos(2), pos(3), pos(4)]);
                        catch
                        end
                    end
                    
                    % Speichern (PNG + FIG)
                    if exist('exportgraphics', 'file')
                        exportgraphics(f_temp, [baseName '.png'], 'Resolution', 150);
                    else
                        saveas(f_temp, [baseName '.png']);
                    end
                    savefig(f_temp, [baseName '.fig']);
                    delete(f_temp);
                    fprintf('OK\n');
                    
                catch innerME
                    fprintf('FEHLER: %s\n', innerME.message);
                    % Wir machen weiter mit dem nächsten
                end
            end
        catch ME
            errordlg(['Fehler beim Batch-Export: ' ME.message], 'Fehler');
        end
        if ishandle(hWait), close(hWait); end
        msgbox(sprintf('Batch Export abgeschlossen.\n%d Dateien verarbeitet.\nOrdner: %s', nExport, outDir), 'Info');
    end

    % --- Heatmap Helper ---
    function updateHeatmapFrame(~, ~)
        % Callback für Slider
        updatePlot();
    end

    function playAnimation(~, ~)
        if isPlaying
            % Pause gedrückt
            isPlaying = false;
            set(hBtnPlay, 'String', 'Play Animation');
            return;
        end

        % Play gedrückt
        isPlaying = true;
        set(hBtnPlay, 'String', 'Pause Animation');

        t_min = get(hSliderTime, 'Min');
        t_max = get(hSliderTime, 'Max');
        step = 0.5; % ms pro Frame
        
        % Start ab aktueller Position (oder von vorne wenn am Ende)
        t_curr = get(hSliderTime, 'Value');
        if t_curr >= t_max - step, t_curr = t_min; end
        
        for t = t_curr:step:t_max
            if ~isPlaying || ~isvalid(f), break; end
            set(hSliderTime, 'Value', t);
            updatePlot();
            drawnow;
            pause(0.05);
        end
        
        % Am Ende zurücksetzen
        if isPlaying && isvalid(hBtnPlay)
            isPlaying = false;
            set(hBtnPlay, 'String', 'Play Animation');
        end
    end

    function dataMap = get_heatmap_data(variante, slot)
        % Lädt alle IRs einer Variante und cached sie
        if slot == 1 && strcmp(heatmapCache.variante1, variante)
            dataMap = heatmapCache.data1; return;
        elseif slot == 2 && strcmp(heatmapCache.variante2, variante)
            dataMap = heatmapCache.data2; return;
        end
        
        % Neu laden: Map mit String Keys ('1'..'15', 'Q1')
        dataMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
        
        % Suche Dateien
        files = dir(fullfile(procDir, sprintf('Proc_%s_*.mat', variante)));
        
        for i = 1:length(files)
            try
                tmp = load(fullfile(files(i).folder, files(i).name), 'Result');
                meta = tmp.Result.meta;
                
                key = '';
                if strcmp(meta.type, 'Source')
                    key = 'Q1';
                else
                    pNum = str2double(meta.position);
                    if ~isnan(pNum)
                        key = num2str(pNum);
                    else
                        digits = regexp(meta.position, '\d+', 'match');
                        if ~isempty(digits), key = digits{1}; end
                    end
                end
                
                if ~isempty(key)
                    dataMap(key) = tmp.Result.time.ir;
                end
            catch
            end
        end
        
        % Cache update
        if slot == 1
            heatmapCache.variante1 = variante; heatmapCache.data1 = dataMap;
        else
            heatmapCache.variante2 = variante; heatmapCache.data2 = dataMap;
        end
    end

    function gridData = calc_heatmap_grid(dataMap, t_ms, fs)
        % Erstellt 4x4 Matrix für imagesc
        gridData = NaN(4, 4);
        win_samp = round(5e-3 * fs); % 10ms Fenster
        idx_start = round(t_ms/1000 * fs) + 1;
        
        % Helper
        get_val = @(k) get_rms_db(dataMap, k, idx_start, win_samp);
        
        % Zeile 1 (Pos 1-4)
        for c=1:4, gridData(1,c) = get_val(num2str(c)); end
        % Zeile 2 (Pos 5-8)
        for c=1:4, gridData(2,c) = get_val(num2str(c+4)); end
        % Zeile 3 (Pos 9-12)
        for c=1:4, gridData(3,c) = get_val(num2str(c+8)); end
        % Zeile 4 (Q1, 13, 14, 15)
        gridData(4,1) = get_val('Q1');
        gridData(4,2) = get_val('13');
        gridData(4,3) = get_val('14');
        gridData(4,4) = get_val('15');
    end
    
    function val = get_rms_db(dataMap, key, idx_start, win_samp)
        val = NaN;
        if isKey(dataMap, key)
            ir = dataMap(key);
            if ~isempty(ir)
                % Indizes berechnen (mit Schutz gegen < 1 für negative Zeiten)
                idx_end = idx_start + win_samp - 1;
                eff_start = max(1, idx_start);
                eff_end = min(length(ir), idx_end);
                
                if eff_start <= eff_end
                    seg = ir(eff_start:eff_end);
                    % RMS über das gesamte Fenster (Nullen annehmen wo keine Daten)
                    rms_val = sqrt(sum(seg.^2) / win_samp);
                    val = 20*log10((rms_val + eps) / FS_global_ref);
                else
                    val = -100;
                end
            else
                val = -100;
            end
        end
    end
    
    function add_heatmap_labels(gridData, threshold)
        layoutLabels = {
            'M1', 'M2', 'M3', 'M4';
            'M5', 'M6', 'M7', 'M8';
            'M9', 'M10', 'M11', 'M12';
            'Q1', 'M13', 'M14', 'M15'
        };
        [rows, cols] = size(layoutLabels);
        for r = 1:rows
            for c = 1:cols
                val = gridData(r,c);
                lbl = layoutLabels{r,c};
                if ~isnan(val) && val > threshold
                    text(c, r, sprintf('%s\n%.0f', lbl, val), ...
                        'HorizontalAlignment', 'center', 'Color', 'k', 'FontSize', 8, 'FontWeight', 'bold');
                else
                    text(c, r, lbl, ...
                        'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5], 'FontSize', 8);
                end
            end
        end
    end

    function s = cleanName(filename)
        % Entfernt 'Proc_' und '.mat' für schönere Legenden
        s = strrep(filename, 'Proc_', '');
        s = strrep(s, 'Time_', '');
        s = strrep(s, '.mat', '');
        s = strrep(s, '_', ' ');
    end

    % --- Helper für Raw Data Loading (falls nicht im Path) ---
    function [S, meta] = load_and_parse_file(filepath)
        [~, fname, ~] = fileparts(filepath);
        S = load(filepath);
        meta = struct();
        meta.filename = fname;
        
        % Regex für Variante und Position
        tokens = regexp(fname, '^(.*?)[_,]Pos[_,]?(\w+)', 'tokens', 'once', 'ignorecase');
        if ~isempty(tokens)
            meta.variante = tokens{1};
            meta.position = tokens{2};
            meta.type = 'Receiver';
        else
            tokens = regexp(fname, '^(.*?)[_,]Quelle', 'tokens', 'once', 'ignorecase');
            if ~isempty(tokens)
                meta.variante = tokens{1};
                meta.position = 'Q1';
                meta.type = 'Source';
            else
                meta.variante = 'Unknown';
                meta.position = '0';
                meta.type = 'Unknown';
            end
        end
    end

    function ir = extract_ir(S)
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
end