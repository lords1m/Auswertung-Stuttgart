function [L_dBFS, L_sum, f_mitten] = calc_terz_spectrum(ir, fs, FS_global)
    % 1. FFT
    N = length(ir);
    N_fft = 2^nextpow2(N); 
    X = fft(ir, N_fft);
    freqs = (0:N_fft-1) * (fs / N_fft);
    
    % Nur positive Frequenzen
    valid_idx = 1:floor(N_fft/2)+1;
    X = X(valid_idx);
    freqs = freqs(valid_idx);
    
    % Energie-Dichte (Parseval-Korrektur)
    X_mag_sq = (abs(X).^2) / N; 
    
    % 2. Terz-Definition (Begrenzt auf 4 kHz - 60 kHz)
    f_mitten = [4000 5000 6300 7100 8000 ...
                9000 10000 11200 12500 14000 16000 18000 20000 22400 25000 ...
                28000 31500 35500 40000 45000 50000 56000];
                
    L_dBFS = NaN(size(f_mitten));
    energy_sum = 0;
    
    for k = 1:length(f_mitten)
        fc = f_mitten(k);
        fl = fc / 2^(1/6);
        fu = fc * 2^(1/6);
        
        if fl > fs/2, break; end
        
        idx = freqs >= fl & freqs <= fu;
        
        if any(idx)
            band_energy = sum(X_mag_sq(idx));
            energy_sum = energy_sum + band_energy;
            L_dBFS(k) = 10 * log10(band_energy / (FS_global^2 + eps));
        else
            L_dBFS(k) = -Inf;
        end
    end
    
    % 3. Summenpegel direkt aus der aufsummierten Energie berechnen
    if energy_sum <= 0
        L_sum = -Inf;
    else
        L_sum = 10 * log10(energy_sum / (FS_global^2 + eps));
    end
end