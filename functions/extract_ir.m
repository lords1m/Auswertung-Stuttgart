function ir = extract_ir(S)
    ir = [];
    field_candidates = {'RiR', 'RIR', 'ir', 'IR', 'aufn', 'audio'};
    
    for i = 1:length(field_candidates)
        fn = field_candidates{i};
        if isfield(S, fn) && isnumeric(S.(fn)) && numel(S.(fn)) > 100
            ir = double(S.(fn)(:));
            return;
        end
    end
    
    fns = fieldnames(S);
    for i = 1:numel(fns)
        if startsWith(fns{i}, '__'), continue; end
        val = S.(fns{i});
        if isnumeric(val) && numel(val) > 1000
            ir = double(val(:));
            return;
        end
    end
end