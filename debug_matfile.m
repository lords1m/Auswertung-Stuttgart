% Debug: Prüfe Inhalt einer .mat-Datei
filename = 'Variante_1_data/Variante_1_neu,Pos_1.mat';
fprintf('Lade: %s\n', filename);
S = load(filename);
fprintf('Felder: %s\n', strjoin(fieldnames(S), ', '));

for f = 1:numel(fieldnames(S))
    fname = fieldnames(S){f};
    v = S.(fname);
    if isnumeric(v)
        fprintf('  %s: numerisch, Size=[%s]\n', fname, num2str(size(v)));
        if numel(v) > 1 && numel(v) < 20
            fprintf('    Werte: %s\n', num2str(v(1:min(5,numel(v)))));
        end
    else
        fprintf('  %s: %s\n', fname, class(v));
    end
end
