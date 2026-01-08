function [ir_trunc, metrics] = truncate_ir(ir)
    % metrics enthält SNR, StartIndex, EndIndex für Nachvollziehbarkeit
    
    % 0. DC-Offset entfernen (wichtig, da konstanter Offset > Schwellwert das Abschneiden verhindert)
    ir = ir - mean(ir);
    
    N = length(ir);
    ir_abs = abs(ir);
    max_amp = max(ir_abs);
    
    if max_amp == 0
        ir_trunc = ir;
        metrics = struct('idx_start',1, 'idx_end',N, 'snr_db',0);
        return;
    end
    
    % 1. Rauschpegel (letzte 10%)
    noise_section = ir_abs(ceil(N*0.9):end);
    if isempty(noise_section), noise_section = ir_abs(end); end % Fallback bei sehr kurzen Signalen
    
    mu_noise = mean(noise_section);
    std_noise = std(noise_section);
    
    % Schwellwert
    % Konservativerer Schwellwert (6 Sigma) und Glättung gegen Spikes
    % Erhöht auf 0.003 (-50dB), um DC-Reste oder Grundrauschen sicher zu unterschreiten
    threshold_end = max(mu_noise + 6*std_noise, max_amp * 0.003); 
    
    % Glättung (Moving Average, ca. 1000 Samples), um einzelne Rauschspitzen zu ignorieren
    window_size = 1000;
    if N > window_size
        env = movmean(ir_abs, window_size);
    else
        env = ir_abs;
    end
    
    % Ende finden: Suche nach dem ersten stabilen Abfall in das Rauschen
    % Statt nach dem allerletzten Punkt über der Schwelle zu suchen (was anfällig für 
    % späte Störgeräusche ist), suchen wir nach dem ersten Zeitraum der Stille nach dem Peak.
    [~, idx_peak] = max(env);
    
    % Wir definieren "Stille" als 20ms (ca. 10k Samples bei 500kHz) unter dem Schwellwert.
    min_silence_samples = 10000; 
    
    if idx_peak < N
        % Suche ab dem Peak vorwärts
        post_peak_env = env(idx_peak:end);
        % Gleitendes Maximum über das Fenster in die Zukunft [0, window]
        local_max_future = movmax(post_peak_env, [0, min_silence_samples]);
        
        % Erster Index, wo das zukünftige Max unter der Schwelle liegt
        idx_silence_start = find(local_max_future < threshold_end, 1, 'first');
        
        if ~isempty(idx_silence_start)
            idx_end = idx_peak + idx_silence_start - 1;
        else
            % Fallback: Kein stabiles Ende gefunden -> nimm das letzte Mal über Schwelle
            idx_end = find(env > threshold_end, 1, 'last');
            if isempty(idx_end), idx_end = N; end
        end
    else
        idx_end = N;
    end
    
    % Puffer hinzufügen (ca. 10ms), um den Ausklang sicher zu haben
    idx_end = min(N, idx_end + 5000);
    
    % Start finden
    threshold_start = max_amp * 0.01; % Startschwelle empfindlicher (1%)
    idx_start = find(ir_abs > threshold_start, 1, 'first');
    
    if isempty(idx_start)
        idx_start = 1; 
    else
        % Pre-Roll: 500 Samples Sicherheit vor dem Impuls
        idx_start = max(1, idx_start - 500);
    end
    
    if idx_end < idx_start, idx_end = N; end
    
    % Zuschneiden
    ir_trunc = ir(idx_start:idx_end);
    
    % Metriken speichern
    metrics.idx_start = idx_start;
    metrics.idx_end = idx_end;
    metrics.original_len = N;
    metrics.snr_db = 20*log10(rms(ir_trunc) / (std_noise + eps));
    metrics.energy = sum(ir_trunc.^2);
    metrics.energy_total = sum(ir.^2);
    metrics.energy_share = metrics.energy / (metrics.energy_total + eps);
end