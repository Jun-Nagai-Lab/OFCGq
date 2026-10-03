function LMM_stats()

    [fileName, filePath] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx,*.xls)'}, ...
        'Select the xlsx file containing LMM/LMM_OWA/LMM_TWA data');
    if isequal(fileName, 0)
        disp('File selection canceled. Exiting.');
        return;
    end
    inputFile = fullfile(filePath, fileName);

    [~, baseName, ~] = fileparts(fileName);
    outputFile = fullfile(filePath, [baseName '_StatsResults.xlsx']);
    if isfile(outputFile)
        delete(outputFile);
    end

    sheets = sheetnames(inputFile);

    for s = 1:numel(sheets)
        sheetName = sheets{s};
        fprintf('--- Processing sheet "%s" ---\n', sheetName);

        raw = readcell(inputFile, 'Sheet', sheetName);
        modeTag  = strtrim(string(raw{1,1}));
        shareTag = strtrim(string(raw{2,1}));

        if ~any(shareTag == ["Share", "no_Share"])
            warning('Sheet "%s": A2 is not "Share"/"no_Share" ("%s"). Defaulting to "no_Share" to be safe.', ...
                sheetName, shareTag);
            shareTag = "no_Share";
        end

        try
            if modeTag == "LMM"
                analyze_LMM_2group(raw, sheetName, shareTag, outputFile);
            elseif modeTag == "LMM_OWA"
                analyze_LMM_OWA(raw, sheetName, shareTag, outputFile);
            elseif modeTag == "LMM_TWA"
                analyze_LMM_TWA(raw, sheetName, shareTag, outputFile);
            else
                warning('Sheet "%s": A1 is not "LMM"/"LMM_OWA"/"LMM_TWA" ("%s"). Skipping.', ...
                    sheetName, modeTag);
                continue;
            end
        catch ME
            warning('An error occurred while processing sheet "%s": %s', sheetName, ME.message);
            if ~isempty(ME.stack)
                fprintf('  -> Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
            end
        end
    end

    fprintf('\nAll sheets processed. Results saved to:\n%s\n', outputFile);
end


function T = parseData(raw)
    header = raw(3, :);
    validCols = ~cellfun(@(x) (ischar(x) && isempty(x)) || ...
        (isstring(x) && x == "") || (isscalar(x) && all(ismissing(x))), header);
    header = header(validCols);
    ncol = numel(header);

    dataRaw = raw(4:end, 1:ncol);
    keepRows = ~cellfun(@(x) (ischar(x) && isempty(x)) || ...
        (isscalar(x) && all(ismissing(x))), dataRaw(:,1));
    dataRaw = dataRaw(keepRows, :);

    varNames = matlab.lang.makeValidName(string(header));
    T = table();
    for c = 1:ncol
        colData = dataRaw(:, c);
        name = varNames(c);
        if any(strcmpi(name, ["Group", "Mouse", "Cell"]))
            T.(name) = string(cellfun(@toStringSafe, colData, 'UniformOutput', false));
        else
            T.(name) = cellfun(@toNumericSafe, colData);
        end
    end
end

function s = toStringSafe(x)
    s = string(x);
end

function v = toNumericSafe(x)
    if isnumeric(x)
        v = double(x);
    elseif ischar(x) || isstring(x)
        v = str2double(x);
    else
        v = NaN;
    end
end

function T = addMouseID(T, shareTag)
    if shareTag == "Share"
        T.MouseID = T.Mouse;
    else
        T.MouseID = T.Group + "_" + T.Mouse;
    end
end

function pAdj = holmBonferroni(p)
    pAdj = nan(size(p));
    validIdx = find(~isnan(p));
    pv = p(validIdx);
    n = numel(pv);
    if n == 0
        return;
    end
    [psorted, idx] = sort(pv);
    padjSorted = nan(n,1);
    for i = 1:n
        padjSorted(i) = min(1, (n - i + 1) * psorted(i));
    end
    for i = 2:n
        padjSorted(i) = max(padjSorted(i), padjSorted(i-1));
    end
    pAdjValid = nan(n,1);
    pAdjValid(idx) = padjSorted;
    pAdj(validIdx) = pAdjValid;
end

function tf = isEssentiallyConstant(vals)
    vals = vals(~isnan(vals));
    if isempty(vals)
        tf = true;
        return;
    end
    scale = max(1, max(abs(vals)));
    tf = (max(vals) - min(vals)) < 1e-10 * scale;
end

function [lme, usedFallback] = fitLMEwithFallback(T, formulaPrimary, formulaFallback)
    usedFallback = false;
    try
        lme = fitlme(T, formulaPrimary);
    catch ME1
        warning('Failed to fit model "%s" (%s). Trying simplified model "%s".', ...
            formulaPrimary, ME1.message, formulaFallback);
        lme = fitlme(T, formulaFallback);
        usedFallback = true;
    end
end

function aov = anovaSatterthwaite(lme)
    aov = anova(lme, 'DFMethod', 'satterthwaite');
end

function nextRow = writeAnovaBlock(outputFile, sheetName, startRow, formulaStr, termLabel, Fval, df1, df2, pval, fallbackNote)
    if nargin < 10
        fallbackNote = '';
    end
    writecell({sprintf('Linear Mixed Model: %s', formulaStr), ''}, ...
        outputFile, 'Sheet', sheetName, 'Range', sprintf('A%d', startRow));

    rows = { ...
        'Effect',               termLabel; ...
        'F',                    Fval; ...
        'df1',                  df1; ...
        'df2 (Satterthwaite)',  df2; ...
        'p',                    pval };
    if ~isempty(fallbackNote)
        rows = [rows; {'Note', fallbackNote}];
    end
    writecell(rows, outputFile, 'Sheet', sheetName, 'Range', sprintf('A%d', startRow + 1));

    nextRow = startRow + 1 + size(rows,1) + 1;
end

function nextRow = writeNotApplicableBlock(outputFile, sheetName, startRow, formulaStr, termLabel, reasonMsg)
    writecell({sprintf('Linear Mixed Model: %s', formulaStr), ''}, ...
        outputFile, 'Sheet', sheetName, 'Range', sprintf('A%d', startRow));
    rows = { ...
        'Effect', termLabel; ...
        'F',   'N/A'; ...
        'df1', 'N/A'; ...
        'df2 (Satterthwaite)', 'N/A'; ...
        'p',   'N/A'; ...
        'Note', reasonMsg };
    writecell(rows, outputFile, 'Sheet', sheetName, 'Range', sprintf('A%d', startRow + 1));
    nextRow = startRow + 1 + size(rows,1) + 1;
end

function [formula, hasCellNesting] = buildLMMFormula(T, valueVar)
    key = string(T.MouseID) + "__" + string(T.Cell);
    [~, ~, ic] = unique(key);
    counts = accumarray(ic, 1);
    hasCellNesting = max(counts) > 1;
    if hasCellNesting
        formula = sprintf('%s ~ Group + (1|MouseID) + (1|MouseID:Cell)', valueVar);
    else
        formula = sprintf('%s ~ Group + (1|MouseID)', valueVar);
    end
end

function analyze_LMM_2group(raw, sheetName, shareTag, outputFile)
    T = parseData(raw);
    T = addMouseID(T, shareTag);

    groups = unique(T.Group, 'stable');
    if numel(groups) ~= 2
        warning('Sheet "%s": Group does not have exactly 2 levels (%d levels found). Using the first 2 in order of appearance.', ...
            sheetName, numel(groups));
        groups = groups(1:2);
    end
    T.Group   = categorical(T.Group, cellstr(groups));
    T.Mouse   = categorical(T.Mouse);
    T.Cell    = categorical(T.Cell);
    T.MouseID = categorical(T.MouseID);

    valueVar = "Value1";
    if ~ismember(valueVar, string(T.Properties.VariableNames))
        numericVars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
        valueVar = string(numericVars{1});
    end

    row = 1;

    if isEssentiallyConstant(T.(valueVar))
        formula = sprintf('%s ~ Group + (1|MouseID)', valueVar);
        row = writeNotApplicableBlock(outputFile, sheetName, row, formula, 'Group', ...
            sprintf('All values of %s are (nearly) identical, so the LMM variance cannot be estimated and the test cannot be performed.', valueVar));
    else
        [formula, hasCellNesting] = buildLMMFormula(T, valueVar);
        if hasCellNesting
            fallbackFormula = sprintf('%s ~ Group + (1|MouseID)', valueVar);
        else
            fallbackFormula = sprintf('%s ~ Group + (1|MouseID:Cell)', valueVar);
        end
        lme = fitLMEwithFallback(T, formula, fallbackFormula);
        aov = anovaSatterthwaite(lme);
        termIdx = find(string(aov.Term) == "Group", 1);
        note = '';
        if hasCellNesting && ~contains(char(lme.Formula), 'MouseID:Cell')
            note = 'A paired structure (nested Cell random effect) was intended but failed to converge; fell back to the simplified model (1|MouseID).';
        end
        row = writeAnovaBlock(outputFile, sheetName, row, char(lme.Formula), 'Group', ...
            aov.FStat(termIdx), aov.DF1(termIdx), aov.DF2(termIdx), aov.pValue(termIdx), note);
    end
    writecell({sprintf('Share flag: %s', shareTag)}, outputFile, 'Sheet', sheetName, ...
        'Range', sprintf('A%d', row));
end

function analyze_LMM_OWA(raw, sheetName, shareTag, outputFile)
    T = parseData(raw);
    T = addMouseID(T, shareTag);

    groups = unique(T.Group, 'stable');
    if numel(groups) < 3
        warning('Sheet "%s": LMM_OWA was specified but Group has only %d levels (3 or more is typically expected). Continuing anyway.', ...
            sheetName, numel(groups));
    end
    T.Group   = categorical(T.Group, cellstr(groups));
    T.Mouse   = categorical(T.Mouse);
    T.Cell    = categorical(T.Cell);
    T.MouseID = categorical(T.MouseID);

    valueVar = "Value1";
    if ~ismember(valueVar, string(T.Properties.VariableNames))
        numericVars = T.Properties.VariableNames(varfun(@isnumeric, T, 'OutputFormat', 'uniform'));
        valueVar = string(numericVars{1});
    end

    row = 1;

    writecell({'=== Omnibus test (One-way LMM) ==='}, outputFile, 'Sheet', sheetName, ...
        'Range', sprintf('A%d', row));
    row = row + 2;

    if isEssentiallyConstant(T.(valueVar))
        formulaOmni = sprintf('%s ~ Group + (1|MouseID)', valueVar);
        row = writeNotApplicableBlock(outputFile, sheetName, row, formulaOmni, 'Group', ...
            sprintf('All values of %s are (nearly) identical, so the LMM variance cannot be estimated and the test cannot be performed.', valueVar));
    else
        [formulaOmni, hasCellNesting] = buildLMMFormula(T, valueVar);
        if hasCellNesting
            fallbackOmni = sprintf('%s ~ Group + (1|MouseID)', valueVar);
        else
            fallbackOmni = sprintf('%s ~ Group + (1|MouseID:Cell)', valueVar);
        end
        lmeOmni = fitLMEwithFallback(T, formulaOmni, fallbackOmni);
        aovOmni = anovaSatterthwaite(lmeOmni);
        termIdx = find(string(aovOmni.Term) == "Group", 1);
        note = '';
        if hasCellNesting && ~contains(char(lmeOmni.Formula), 'MouseID:Cell')
            note = 'A paired structure (nested Cell random effect) was intended but failed to converge; fell back to the simplified model (1|MouseID).';
        end
        row = writeAnovaBlock(outputFile, sheetName, row, char(lmeOmni.Formula), 'Group', ...
            aovOmni.FStat(termIdx), aovOmni.DF1(termIdx), aovOmni.DF2(termIdx), aovOmni.pValue(termIdx), note);
    end

    writecell({sprintf('Share flag: %s', shareTag)}, outputFile, 'Sheet', sheetName, ...
        'Range', sprintf('A%d', row));
    row = row + 2;

    writecell({'=== Post-hoc comparisons (pairwise Group comparisons) ==='}, ...
        outputFile, 'Sheet', sheetName, 'Range', sprintf('A%d', row));
    row = row + 2;

    nGroups = numel(groups);
    pairIdx = nchoosek(1:nGroups, 2);
    nPairs = size(pairIdx, 1);

    Group1_ = strings(nPairs,1); Group2_ = strings(nPairs,1);
    FvalArr = nan(nPairs,1); df1Arr = nan(nPairs,1); df2Arr = nan(nPairs,1);
    pRaw = nan(nPairs,1);

    for k = 1:nPairs
        gA = groups(pairIdx(k,1));
        gB = groups(pairIdx(k,2));
        sub = T(T.Group == gA | T.Group == gB, :);
        sub.Group = categorical(string(sub.Group), cellstr([gA; gB]));

        formulaPair = sprintf('%s ~ Group + (1|MouseID)', valueVar);

        writecell({sprintf('--- %s vs %s ---', gA, gB)}, outputFile, 'Sheet', sheetName, ...
            'Range', sprintf('A%d', row));
        row = row + 1;

        if isEssentiallyConstant(sub.(valueVar))
            Group1_(k) = gA; Group2_(k) = gB;
            FvalArr(k) = NaN; df1Arr(k) = NaN; df2Arr(k) = NaN; pRaw(k) = NaN;
            row = writeNotApplicableBlock(outputFile, sheetName, row, formulaPair, 'Group', ...
                sprintf('All values of %s are (nearly) identical in both groups, so the LMM variance cannot be estimated and the test cannot be performed.', valueVar));
        else
            [formulaPair, hasCellNesting] = buildLMMFormula(sub, valueVar);
            if hasCellNesting
                fallbackPair = sprintf('%s ~ Group + (1|MouseID)', valueVar);
            else
                fallbackPair = sprintf('%s ~ Group + (1|MouseID:Cell)', valueVar);
            end
            lmePair = fitLMEwithFallback(sub, formulaPair, fallbackPair);
            aovPair = anovaSatterthwaite(lmePair);
            tIdx = find(string(aovPair.Term) == "Group", 1);

            Group1_(k) = gA; Group2_(k) = gB;
            FvalArr(k) = aovPair.FStat(tIdx);
            df1Arr(k)  = aovPair.DF1(tIdx);
            df2Arr(k)  = aovPair.DF2(tIdx);
            pRaw(k)    = aovPair.pValue(tIdx);

            note = '';
            if hasCellNesting && ~contains(char(lmePair.Formula), 'MouseID:Cell')
                note = 'A paired structure (nested Cell random effect) was intended but failed to converge; fell back to the simplified model (1|MouseID).';
            end
            row = writeAnovaBlock(outputFile, sheetName, row, char(lmePair.Formula), 'Group', ...
                FvalArr(k), df1Arr(k), df2Arr(k), pRaw(k), note);
        end
    end

    pAdj = holmBonferroni(pRaw);
    summaryTable = table(Group1_, Group2_, FvalArr, df1Arr, df2Arr, pRaw, pAdj, ...
        'VariableNames', {'Group1','Group2','F','df1','df2_Satterthwaite','pValue_raw','pValue_HolmAdj'});

    writecell({'=== Post-hoc summary (Holm-Bonferroni corrected) ==='}, ...
        outputFile, 'Sheet', sheetName, 'Range', sprintf('A%d', row));
    row = row + 1;
    writetable(summaryTable, outputFile, 'Sheet', sheetName, 'Range', sprintf('A%d', row));
end

function analyze_LMM_TWA(raw, sheetName, shareTag, outputFile)
    T = parseData(raw);
    T = addMouseID(T, shareTag);

    groups = unique(T.Group, 'stable');
    T.Group   = categorical(T.Group, cellstr(groups));
    T.Mouse   = categorical(T.Mouse);
    T.Cell    = categorical(T.Cell);
    T.MouseID = categorical(T.MouseID);

    allVars = string(T.Properties.VariableNames);
    valueVars = allVars(startsWith(allVars, "Value"));
    if isempty(valueVars)
        error('Sheet "%s": No Value* columns found for LMM_TWA analysis.', sheetName);
    end
    nVal = numel(valueVars);
    nRows = height(T);

    Group_l   = repmat(T.Group,   nVal, 1);
    Mouse_l   = repmat(T.Mouse,   nVal, 1);
    Cell_l    = repmat(T.Cell,    nVal, 1);
    MouseID_l = repmat(T.MouseID, nVal, 1);
    ValueName_l = strings(nRows*nVal, 1);
    Value_l = nan(nRows*nVal, 1);
    for v = 1:nVal
        rIdx = (v-1)*nRows+1 : v*nRows;
        ValueName_l(rIdx) = valueVars(v);
        Value_l(rIdx) = T.(valueVars(v));
    end
    Long = table(Group_l, Mouse_l, Cell_l, MouseID_l, ...
        categorical(ValueName_l, cellstr(valueVars)), Value_l, ...
        'VariableNames', {'Group','Mouse','Cell','MouseID','ValueName','Value'});

    row = 1;

    formulaOmni = 'Value ~ Group*ValueName + (1|MouseID) + (1|MouseID:Cell)';
    fallbackOmni = 'Value ~ Group*ValueName + (1|MouseID:Cell)';
    lmeOmni = fitLMEwithFallback(Long, formulaOmni, fallbackOmni);
    aovOmni = anovaSatterthwaite(lmeOmni);

    writecell({'=== Omnibus test (Two-way LMM) ==='}, outputFile, 'Sheet', sheetName, ...
        'Range', sprintf('A%d', row));
    row = row + 2;
    for i = 1:height(aovOmni)
        termLabel = char(string(aovOmni.Term(i)));
        row = writeAnovaBlock(outputFile, sheetName, row, char(lmeOmni.Formula), termLabel, ...
            aovOmni.FStat(i), aovOmni.DF1(i), aovOmni.DF2(i), aovOmni.pValue(i));
    end

    writecell({sprintf('Share flag: %s', shareTag)}, outputFile, 'Sheet', sheetName, ...
        'Range', sprintf('A%d', row));
    row = row + 2;

    writecell({'=== Post-hoc comparisons (Group vs Group at each Value) ==='}, ...
        outputFile, 'Sheet', sheetName, 'Range', sprintf('A%d', row));
    row = row + 2;

    pRaw = nan(nVal, 1);
    FvalArr = nan(nVal,1); df1Arr = nan(nVal,1); df2Arr = nan(nVal,1);

    for v = 1:nVal
        sub = Long(Long.ValueName == valueVars(v), :);
        sub.Group = categorical(string(sub.Group), cellstr(groups));

        writecell({sprintf('--- Value: %s ---', valueVars(v))}, outputFile, 'Sheet', sheetName, ...
            'Range', sprintf('A%d', row));
        row = row + 1;

        if isEssentiallyConstant(sub.Value)
            FvalArr(v) = NaN; df1Arr(v) = NaN; df2Arr(v) = NaN; pRaw(v) = NaN;
            formulaSub = 'Value ~ Group + (1|MouseID)';
            row = writeNotApplicableBlock(outputFile, sheetName, row, formulaSub, 'Group', ...
                sprintf('All values of %s are (nearly) identical across all rows (e.g., a normalization reference time point), so the LMM variance cannot be estimated and the test cannot be performed.', valueVars(v)));
        else
            [formulaSub, hasCellNesting] = buildLMMFormula(sub, 'Value');
            if hasCellNesting
                fallbackSub = 'Value ~ Group + (1|MouseID)';
            else
                fallbackSub = 'Value ~ Group + (1|MouseID:Cell)';
            end
            lmeSub = fitLMEwithFallback(sub, formulaSub, fallbackSub);
            aovSub = anovaSatterthwaite(lmeSub);
            termIdx = find(string(aovSub.Term) == "Group", 1);

            FvalArr(v) = aovSub.FStat(termIdx);
            df1Arr(v)  = aovSub.DF1(termIdx);
            df2Arr(v)  = aovSub.DF2(termIdx);
            pRaw(v)    = aovSub.pValue(termIdx);

            note = '';
            if hasCellNesting && ~contains(char(lmeSub.Formula), 'MouseID:Cell')
                note = 'A paired structure (nested Cell random effect) was intended but failed to converge; fell back to the simplified model (1|MouseID).';
            end
            row = writeAnovaBlock(outputFile, sheetName, row, char(lmeSub.Formula), 'Group', ...
                FvalArr(v), df1Arr(v), df2Arr(v), pRaw(v), note);
        end
    end

    pAdj = holmBonferroni(pRaw);
    summaryTable = table(valueVars(:), FvalArr, df1Arr, df2Arr, pRaw, pAdj, ...
        'VariableNames', {'Value','F','df1','df2_Satterthwaite','pValue_raw','pValue_HolmAdj'});

    writecell({'=== Post-hoc summary (Holm-Bonferroni corrected) ==='}, ...
        outputFile, 'Sheet', sheetName, 'Range', sprintf('A%d', row));
    row = row + 1;
    writetable(summaryTable, outputFile, 'Sheet', sheetName, 'Range', sprintf('A%d', row));
end