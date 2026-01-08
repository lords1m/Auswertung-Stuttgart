function [S, meta] = load_and_parse_file(filepath)
% LOAD_AND_PARSE_FILE Lädt .mat Datei und extrahiert Metadaten
% UPDATE: "Sanitizing" für macOS -> Ersetzt Kommata/Leerzeichen durch _

    if ~exist(filepath, 'file')
        error('Datei nicht gefunden: %s', filepath);
    end

    % 1. Datei laden
    try
        S = load(filepath);
    catch ME
        warning('Konnte Datei nicht laden: %s', filepath);
        S = []; meta = []; return;
    end

    % 2. Dateinamen analysieren (ohne Pfad und Extension)
    [~, fname, ~] = fileparts(filepath);
    meta.filename = fname;
    meta.variante = 'Unknown';
    meta.position = '';
    meta.type = 'Unknown';

    % --- Strategie A: Position finden (Suche nach "Pos" am Ende) ---
    % Regex sucht von RECHTS: Trenner + "Pos" + Trenner + ID + Ende($)
    [tokens, startIdx] = regexp(fname, '(?i)[_,;\. ]+pos[_\- ]*([A-Za-z0-9_]+)$', 'tokens', 'start', 'once');
    
    if ~isempty(tokens)
        rawPos = tokens{1};
        rawVariante = fname(1:startIdx-1);
        
        meta.position = sanitize_string(rawPos);
        meta.variante = sanitize_string(rawVariante);
        meta.type = 'Mic';
        return;
    end
    
    % --- Strategie B: Quelle finden ---
    [tokens_src, startIdx_src] = regexp(fname, '(?i)[_,;\. ]+quelle.*$', 'tokens', 'start', 'once');
    
    if ~isempty(startIdx_src)
        rawVariante = fname(1:startIdx_src-1);
        
        meta.position = '0';
        meta.variante = sanitize_string(rawVariante);
        meta.type = 'Source';
        return;
    end
    
    warning('Dateiname entspricht keinem bekannten Muster: %s', fname);
end

function cleanStr = sanitize_string(str)
    % Hilfsfunktion: Macht Strings sicher für Dateinamen
    % 1. Ersetze Kommas, Punkte, Leerzeichen, Bindestriche durch Unterstrich
    cleanStr = regexprep(str, '[,;\. \-]', '_');
    
    % 2. Entferne doppelte Unterstriche (z.B. "__" -> "_")
    cleanStr = regexprep(cleanStr, '_+', '_');
    
    % 3. Entferne Unterstriche am Anfang oder Ende
    cleanStr = regexprep(cleanStr, '^_+|_+$', '');
    
    % 4. Leerzeichen sicherheitshalber entfernen (falls Regex was übersehen hat)
    cleanStr = strtrim(cleanStr);
end