function PeriBlink_Calcium()

params.caTimeUnit      = 'sec';
params.fps             = 30;
params.upperY_col      = 15;
params.lowerY_col      = 18;
params.win_sec_dee     = 30;

params.f0Window        = [5 25] * 60;
params.f0Percentile    = 25;
params.analysisWindow  = [5 25] * 60;

params.blinkDEEThreshold      = 0.5;
params.blinkDurationThreshold = 0;

params.periBlinkWindowSec     = 10;
params.periBlinkZeroAtOnset   = false;
params.timeCourseBinWidthSec  = 1;

params.smoothWinSec_dee       = 0;
params.smoothWinSec_caDFF     = 0;

[fosPosFile, fosPosPath] = uigetfile('*.xlsx', 'Select the Ca signal file for Fos+ cells');
if isequal(fosPosFile, 0), disp('Cancelled'); return; end
[fosNegFile, fosNegPath] = uigetfile('*.xlsx', 'Select the Ca signal file for Fos- cells');
if isequal(fosNegFile, 0), disp('Cancelled'); return; end
dlcFolder = uigetdir(pwd, 'Select the folder containing the DLC output csv files');
if isequal(dlcFolder, 0), disp('Cancelled'); return; end

outputFolder = fullfile(fosPosPath, 'PeriBlink_Output');
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end
tiffPath  = fullfile(outputFolder, 'PeriBlink_Calcium_Traces.tif');
excelPath = fullfile(outputFolder, 'Analysis_Data.xlsx');
if exist(tiffPath, 'file'),  delete(tiffPath);  end
if exist(excelPath, 'file'), delete(excelPath); end
fprintf('Output folder: %s\n\n', outputFolder);

meta = table('Size', [0 4], 'VariableTypes', {'string', 'string', 'double', 'string'}, ...
    'VariableNames', {'Analysis', 'Metric', 'Value', 'Details'});
meta = appendParamsToMeta(meta, params);
meta = appendMethodToMeta(meta);

fprintf('[Step 1] Loading calcium data...\n');
fosPos = loadCaFile(fullfile(fosPosPath, fosPosFile), params.caTimeUnit);
fosNeg = loadCaFile(fullfile(fosNegPath, fosNegFile), params.caTimeUnit);
fprintf('  Fos+: %d cells / Fos-: %d cells\n\n', size(fosPos.Ca, 1), size(fosNeg.Ca, 1));

meta = appendMetaRow(meta, 'Sample_Size', 'N_FosPos_Cells', size(fosPos.Ca, 1), '');
meta = appendMetaRow(meta, 'Sample_Size', 'N_FosNeg_Cells', size(fosNeg.Ca, 1), '');
meta = appendMetaRow(meta, 'Sample_Size', 'N_Mice', numel(unique([fosPos.mouseLabels; fosNeg.mouseLabels])), '');

fprintf('[Step 2] Computing DFF...\n');
[tPos, dffPos] = computeDFF(fosPos, params);
[tNeg, dffNeg] = computeDFF(fosNeg, params);

fprintf('[Step 3] Detecting blinks based on DEE_Inverted...\n');
deeMap      = computeDEEInvertedMap(dlcFolder, params);
deeMap      = restrictMapToWindow(deeMap, params.analysisWindow);
blinkEvents = detectBlinkEvents(deeMap, params.blinkDEEThreshold, params.fps);
blinkEvents.IncludedInPeriAnalysis = blinkEvents.Duration_sec >= params.blinkDurationThreshold;
longBlinks  = blinkEvents(blinkEvents.IncludedInPeriAnalysis, :);
fprintf('  -> Total blinks detected: %d / Duration >= %.2f sec: %d\n', ...
    height(blinkEvents), params.blinkDurationThreshold, height(longBlinks));

meta = appendMetaRow(meta, 'Blink_Detection', 'N_Total_Blinks', height(blinkEvents), ...
    sprintf('DEE_Inverted >= %.3f; analysis window = %g-%g sec', ...
    params.blinkDEEThreshold, params.analysisWindow(1), params.analysisWindow(2)));
meta = appendMetaRow(meta, 'Blink_Detection', 'N_Long_Blinks', height(longBlinks), ...
    sprintf('Duration >= %.3f sec (used for Peri-Blink analysis)', params.blinkDurationThreshold));

timeCourseRaw = table();
if isempty(longBlinks)
    warning('No blink events met the criteria; skipping Peri-Blink analysis.');
else
    fprintf('[Step 4] Extracting Peri-Blink traces and analyzing time course...\n');

    [traceFosPos, relGrid] = extractPeriBlinkTraces(tPos, dffPos, fosPos.mouseLabels, longBlinks, ...
        params.periBlinkWindowSec, params.periBlinkZeroAtOnset);
    traceFosNeg            = extractPeriBlinkTraces(tNeg, dffNeg, fosNeg.mouseLabels, longBlinks, ...
        params.periBlinkWindowSec, params.periBlinkZeroAtOnset);

    binWidth = params.timeCourseBinWidthSec;
    [binCenters, binPos] = binTraces(relGrid, traceFosPos, binWidth);
    [~,          binNeg] = binTraces(relGrid, traceFosNeg, binWidth);
    nBins = numel(binCenters);

    okPos = ~all(isnan(traceFosPos), 2);
    okNeg = ~all(isnan(traceFosNeg), 2);
    mousePos = fosPos.mouseLabels(okPos); cellPos = fosPos.cellLabels(okPos); binPos = binPos(okPos, :);
    mouseNeg = fosNeg.mouseLabels(okNeg); cellNeg = fosNeg.cellLabels(okNeg); binNeg = binNeg(okNeg, :);

    colNames = cell(1, nBins);
    for b = 1:nBins
        if binCenters(b) < 0
            colNames{b} = matlab.lang.makeValidName(sprintf('t_neg%gsec', abs(binCenters(b))));
        else
            colNames{b} = matlab.lang.makeValidName(sprintf('t_pos%gsec', binCenters(b)));
        end
    end
    idTable = table([repmat("Fos+", numel(mousePos), 1); repmat("Fos-", numel(mouseNeg), 1)], ...
        [string(mousePos(:)); string(mouseNeg(:))], ...
        [string(cellPos(:));  string(cellNeg(:))], ...
        'VariableNames', {'Group', 'MouseID', 'CellID'});
    timeCourseRaw = [idTable, array2table([binPos; binNeg], 'VariableNames', colNames)];

    [tcValue, tcGroup, tcMouse, tcCell, tcBin] = toLongFormat(binPos, mousePos, cellPos, binNeg, mouseNeg, cellNeg);

    pGroup = NaN; pTime = NaN; pInter = NaN;
    if ~isempty(tcValue) && numel(unique(tcGroup)) == 2
        [pGroup, pTime, pInter] = runTwoWayLMM(tcValue, tcGroup, tcMouse, tcCell, tcBin);
    end
    fprintf('  Two-way LMM: Group p = %.4g, Time p = %.4g, Group x Time p = %.4g\n', pGroup, pTime, pInter);
    lmmFormula = sprintf('Value ~ Group*TimeBin + (1|Mouse) + (1|Mouse:Cell); %g-sec bins', binWidth);
    meta = appendMetaRow(meta, 'PeriBlink_TimeCourse_TwoWayLMM', 'Group_p',        pGroup, lmmFormula);
    meta = appendMetaRow(meta, 'PeriBlink_TimeCourse_TwoWayLMM', 'Time_p',         pTime,  lmmFormula);
    meta = appendMetaRow(meta, 'PeriBlink_TimeCourse_TwoWayLMM', 'Group_x_Time_p', pInter, lmmFormula);

    pRaw = nan(nBins, 1);
    for b = 1:nBins
        idx = (tcBin == b);
        if sum(tcGroup(idx) == "Fos+") >= 2 && sum(tcGroup(idx) == "Fos-") >= 2
            pRaw(b) = runGroupLMM(tcValue(idx), tcGroup(idx), tcMouse(idx));
        end
    end
    pHolm = holmBonferroni(pRaw);
    for b = 1:nBins
        name = sprintf('PeriBlink_PostHoc_t%+gsec', binCenters(b));
        meta = appendMetaRow(meta, name, 'Group_p_LMM_raw', pRaw(b), ...
            sprintf('Value ~ Group + (1|Mouse), within %g-sec bin (uncorrected)', binWidth));
        meta = appendMetaRow(meta, name, 'Group_p_LMM_Holm', pHolm(b), ...
            sprintf('Holm-Bonferroni corrected across %d bins', nBins));
    end

    exportMousePages(relGrid, traceFosPos, fosPos.mouseLabels, traceFosNeg, fosNeg.mouseLabels, tiffPath);

    [mPos, sPos] = meanSEM(binPos);
    [mNeg, sNeg] = meanSEM(binNeg);
    appendTimeCoursePage(binCenters, mPos, sPos, mNeg, sNeg, pHolm, pGroup, pInter, tiffPath, binWidth, ...
        sprintf('All Cells, n=%d/%d', size(binPos, 1), size(binNeg, 1)));

    [mPosM, sPosM, nMicePos] = meanSEMByMouse(binPos, mousePos);
    [mNegM, sNegM, nMiceNeg] = meanSEMByMouse(binNeg, mouseNeg);
    appendTimeCoursePage(binCenters, mPosM, sPosM, mNegM, sNegM, pHolm, pGroup, pInter, tiffPath, binWidth, ...
        sprintf('Mouse-level Mean \\pm SEM, n=%d/%d mice', nMicePos, nMiceNeg));
end

fprintf('[Step 5] Writing Excel file...\n');
writetable(blinkEvents, excelPath, 'Sheet', 'Blink_Events');
if ~isempty(timeCourseRaw)
    writetable(timeCourseRaw, excelPath, 'Sheet', 'PeriBlink_TimeCourse_Data');
end
writetable(meta, excelPath, 'Sheet', 'Meta data');

fprintf('\n%s\n[Done] Output folder: %s\n%s\n', repmat('=', 1, 60), outputFolder, repmat('=', 1, 60));
end

function S = loadCaFile(filename, timeUnit)
raw = readcell(filename);

timeRaw = raw(3:end, 1);
validRow = true(numel(timeRaw), 1);
for i = 1:numel(timeRaw)
    v = timeRaw{i};
    if isempty(v) || (isscalar(v) && ismissing(v)), validRow(i) = false;
    elseif (ischar(v) || isstring(v)) && strlength(strtrim(string(v))) == 0, validRow(i) = false;
    end
end
timeVec = cellToNumeric(timeRaw(validRow), 3, 1, filename);
if strcmpi(timeUnit, 'min'), timeVec = timeVec * 60; end

mouseRaw = string(raw(1, 2:end));
cellRaw  = string(raw(2, 2:end));
validCol = ~ismissing(mouseRaw) & ~ismissing(cellRaw) & (strtrim(mouseRaw) ~= "") & (strtrim(cellRaw) ~= "");
mouseLabels = strtrim(mouseRaw(validCol));
cellLabels  = strtrim(cellRaw(validCol));
for i = 1:numel(cellLabels)
    tok = regexp(cellLabels(i), "['""](.*?)['""]", 'tokens');
    if ~isempty(tok), cellLabels(i) = string(tok{1}{1}); end
end

block = raw(3:end, 2:end);
block = block(validRow, validCol);

S.time        = timeVec(:)';
S.Ca          = cellToNumeric(block, 3, 2, filename).';
S.mouseLabels = mouseLabels(:);
S.cellLabels  = cellLabels(:);
end

function M = cellToNumeric(block, rowStart, colStart, context)
[nr, nc] = size(block);
M = nan(nr, nc);
bad = {};
for i = 1:nr
    for j = 1:nc
        v = block{i, j};
        if isnumeric(v) && isscalar(v)
            M(i, j) = v;
        elseif islogical(v)
            M(i, j) = double(v);
        elseif isempty(v)
            M(i, j) = NaN;
        elseif ischar(v) || isstring(v)
            s = strtrim(string(v));
            if strlength(s) == 0
                M(i, j) = NaN;
            else
                n = str2double(s);
                if isnan(n), bad{end+1} = sprintf('  Row %d, Col %d = "%s"', i+rowStart-1, j+colStart-1, s);
                else, M(i, j) = n;
                end
            end
        else
            bad{end+1} = sprintf('  Row %d, Col %d : class "%s"', i+rowStart-1, j+colStart-1, class(v));
        end
    end
end
if ~isempty(bad), error('Numeric conversion failed in %s:\n%s', context, strjoin(bad, '\n')); end
end

function [tWin, dff] = computeDFF(S, params)
inF0 = S.time >= params.f0Window(1) & S.time <= params.f0Window(2);
CaF0 = S.Ca(:, inF0);
F0 = nan(size(CaF0, 1), 1);
for c = 1:size(CaF0, 1)
    v = CaF0(c, ~isnan(CaF0(c, :)));
    if ~isempty(v), F0(c) = prctile(v, params.f0Percentile); end
end

inWin = S.time >= params.analysisWindow(1) & S.time <= params.analysisWindow(2);
tWin  = S.time(inWin);
dff   = (S.Ca(:, inWin) - F0) ./ F0;
dff   = movingAverageSmooth(dff, mean(diff(tWin)), params.smoothWinSec_caDFF);
end

function y = movingAverageSmooth(x, dt, winSec)
if isempty(winSec) || winSec <= 0, y = x; return; end
n = max(1, round(winSec / dt));
if isvector(x)
    y = movmean(x(:), n, 'omitnan');
    if isrow(x), y = y'; end
else
    y = movmean(x, n, 2, 'omitnan');
end
end

function deeMap = computeDEEInvertedMap(dlcFolder, params)
csvList = dir(fullfile(dlcFolder, '*.csv'));
deeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = 1:numel(csvList)
    [~, stem] = fileparts(csvList(k).name);
    parts = strsplit(stem, '_');
    mouseID = parts{1};
    try
        data = readmatrix(fullfile(dlcFolder, csvList(k).name), 'NumHeaderLines', 3);
        eyeOpening = data(:, params.lowerY_col) - data(:, params.upperY_col);
        eyeOpening(isnan(eyeOpening)) = 0;
        dee = -calculateDEE(eyeOpening, params.fps, params.win_sec_dee);
        dee = movingAverageSmooth(dee, 1/params.fps, params.smoothWinSec_dee);
        deeMap(mouseID) = struct('time', ((0:numel(dee)-1) / params.fps)', 'value', dee);
    catch ME
        fprintf('  -> [ERROR] Failed to process %s: %s\n', csvList(k).name, ME.message);
    end
end
end

function dee = calculateDEE(signal, fps, winSec)
signal = signal(:);
signal(isnan(signal)) = 0;
movMed = movmedian(signal, winSec * fps * 2 + 1, 'omitnan');
dee = (signal - movMed) ./ movMed;
end

function mapOut = restrictMapToWindow(mapIn, tRange)
mapOut = containers.Map('KeyType', 'char', 'ValueType', 'any');
for id = keys(mapIn)
    s = mapIn(id{1});
    idx = s.time >= tRange(1) & s.time <= tRange(2);
    mapOut(id{1}) = struct('time', s.time(idx), 'value', s.value(idx));
end
end

function T = detectBlinkEvents(deeMap, threshold, fps)
dt = 1 / fps;
mouse = strings(0, 1); tStart = zeros(0, 1); tEnd = zeros(0, 1);
for id = keys(deeMap)
    s = deeMap(id{1});
    d = diff([false; s.value(:) > threshold; false]);
    iStart = find(d == 1);
    iEnd   = find(d == -1) - 1;
    for e = 1:numel(iStart)
        mouse(end+1, 1)  = string(id{1});
        tStart(end+1, 1) = s.time(iStart(e));
        tEnd(end+1, 1)   = s.time(iEnd(e)) + dt;
    end
end
T = table(mouse, tStart, tEnd, tEnd - tStart, ...
    'VariableNames', {'MouseID', 'StartTime_sec', 'EndTime_sec', 'Duration_sec'});
if height(T) > 0, T = sortrows(T, {'MouseID', 'StartTime_sec'}); end
end

function [periTrace, relGrid] = extractPeriBlinkTraces(timeVec, dff, mouseLabels, blinkEvents, winSec, zeroAtOnset)
dt = mean(diff(timeVec));
relGrid = (-winSec:dt:winSec)';
[~, onsetIdx] = min(abs(relGrid));
nCell = size(dff, 1);
periTrace = nan(nCell, numel(relGrid));

for c = 1:nCell
    starts = blinkEvents.StartTime_sec(blinkEvents.MouseID == string(mouseLabels(c)));
    segs = [];
    for tStart = starts(:)'
        tQuery = tStart + relGrid;
        if tQuery(1) < timeVec(1) || tQuery(end) > timeVec(end), continue; end
        seg = interp1(timeVec, dff(c, :), tQuery, 'linear');
        if zeroAtOnset, seg = seg - seg(onsetIdx); end
        segs = [segs; seg(:)'];
    end
    if ~isempty(segs), periTrace(c, :) = mean(segs, 1, 'omitnan'); end
end
end

function [binCenters, binned] = binTraces(relGrid, periTrace, binWidth)
winSec = floor(max(abs(relGrid)));
binCenters = (-winSec:binWidth:winSec)';
binned = nan(size(periTrace, 1), numel(binCenters));
for b = 1:numel(binCenters)
    inBin = relGrid >= (binCenters(b) - binWidth/2) & relGrid < (binCenters(b) + binWidth/2 + 1e-9);
    if any(inBin), binned(:, b) = mean(periTrace(:, inBin), 2, 'omitnan'); end
end
end

function [value, group, mouse, cell, binID] = toLongFormat(binPos, mousePos, cellPos, binNeg, mouseNeg, cellNeg)
value = []; group = strings(0, 1); mouse = strings(0, 1); cell = strings(0, 1); binID = [];
groups = {"Fos+", binPos, mousePos, cellPos; "Fos-", binNeg, mouseNeg, cellNeg};
for b = 1:size(binPos, 2)
    for g = 1:2
        v  = groups{g, 2}(:, b);
        ok = ~isnan(v);
        n  = sum(ok);
        value = [value; v(ok)];
        group = [group; repmat(groups{g, 1}, n, 1)];
        mouse = [mouse; string(groups{g, 3}(ok))];
        cell  = [cell;  string(groups{g, 4}(ok))];
        binID = [binID; repmat(b, n, 1)];
    end
end
end

function [m, s] = meanSEM(X)
m = nan(size(X, 2), 1); s = nan(size(X, 2), 1);
for b = 1:size(X, 2)
    v = X(~isnan(X(:, b)), b);
    if ~isempty(v), m(b) = mean(v); s(b) = std(v) / sqrt(numel(v)); end
end
end

function [m, s, nMice] = meanSEMByMouse(X, mouseLabels)
u = unique(string(mouseLabels));
nMice = numel(u);
M = nan(nMice, size(X, 2));
for i = 1:nMice
    M(i, :) = mean(X(string(mouseLabels) == u(i), :), 1, 'omitnan');
end
[m, s] = meanSEM(M);
end

function p = runGroupLMM(value, group, mouse)
T = table(value(:), categorical(string(group(:))), categorical(string(mouse(:))), ...
    'VariableNames', {'Value', 'Group', 'Mouse'});
mdl = fitlme(T, 'Value ~ Group + (1|Mouse)');
p = mdl.Coefficients.pValue(2);
end

function [pGroup, pTime, pInter] = runTwoWayLMM(value, group, mouse, cell, binID)
T = table(value(:), categorical(string(group(:))), categorical(string(mouse(:))), ...
    categorical(string(cell(:))), categorical(string(binID(:))), ...
    'VariableNames', {'Value', 'Group', 'Mouse', 'Cell', 'TimeBin'});
mdl = fitlme(T, 'Value ~ Group*TimeBin + (1|Mouse) + (1|Mouse:Cell)');
a = anova(mdl);
pGroup = a.pValue(strcmp(a.Term, 'Group'));
pTime  = a.pValue(strcmp(a.Term, 'TimeBin'));
pInter = a.pValue(strcmp(a.Term, 'Group:TimeBin'));
end

function pAdj = holmBonferroni(p)
pAdj = nan(size(p));
idx = find(~isnan(p(:)));
if isempty(idx), return; end
m = numel(idx);
[ps, order] = sort(p(idx));
adj = min(cummax(ps(:) .* (m - (1:m)' + 1)), 1);
tmp = nan(m, 1);
tmp(order) = adj;
pAdj(idx) = tmp;
end

function appendImage(fig, savePath)
drawnow;
im = frame2im(getframe(fig));
if exist(savePath, 'file'), mode = 'append'; else, mode = 'overwrite'; end
imwrite(im, savePath, 'TIFF', 'WriteMode', mode);
end

function [mouseMean, uMouse] = meanByMouse(traces, mouseLabels)
uMouse = unique(string(mouseLabels), 'stable');
mouseMean = nan(size(traces, 2), numel(uMouse));
for i = 1:numel(uMouse)
    mouseMean(:, i) = mean(traces(string(mouseLabels) == uMouse(i), :), 1, 'omitnan')';
end
end

function plotTraceAxes(ax, relGrid, y, color, lineWidth, titleStr, titleSize, titleWeight)
hold(ax, 'on');
if isempty(y)
    text(ax, 0.5, 0.5, 'No Data', 'Units', 'normalized', 'HorizontalAlignment', 'center');
else
    plot(ax, relGrid, y, 'Color', color, 'LineWidth', lineWidth);
end
yl = ylim(ax); yl = [min(yl(1), 0), max(yl(2), 0)];
ylim(ax, yl);
line(ax, [0 0], yl, 'Color', [0.3 0.3 0.3], 'LineStyle', '--');
xlim(ax, [min(relGrid) max(relGrid)]);
ylabel(ax, 'DFF'); grid(ax, 'on');
title(ax, titleStr, 'FontName', 'Arial', 'FontSize', titleSize, 'FontWeight', titleWeight, 'Interpreter', 'none');
set(ax, 'FontName', 'Arial', 'FontSize', 8, 'Box', 'off', 'TickDir', 'out');
end

function exportMousePages(relGrid, tracePos, mousePos, traceNeg, mouseNeg, savePath)
colPos = [0.85 0.1 0.1]; colNeg = [0.1 0.1 0.85];
[mmPos, uPos] = meanByMouse(tracePos, mousePos);
[mmNeg, uNeg] = meanByMouse(traceNeg, mouseNeg);
okPos = ~all(isnan(tracePos), 2);
okNeg = ~all(isnan(traceNeg), 2);
allMice = unique([string(mousePos); string(mouseNeg)], 'stable');

fig = figure('Units', 'pixels', 'Position', [100 100 700 500], 'Color', 'w', 'Visible', 'off');
cleanup = onCleanup(@() close(fig));

for m = 1:numel(allMice)
    id = allMice(m);
    clf(fig);
    iP = find(uPos == id, 1); iN = find(uNeg == id, 1);
    yP = []; if ~isempty(iP), yP = mmPos(:, iP); end
    yN = []; if ~isempty(iN), yN = mmNeg(:, iN); end

    ax1 = subplot(2, 1, 1, 'Parent', fig);
    plotTraceAxes(ax1, relGrid, yP, colPos, 1.5, ...
        sprintf('Mouse: %s | Fos+ Peri-Blink Mean Ca (n=%d cells)', id, sum(string(mousePos) == id & okPos)), 9, 'normal');
    ax2 = subplot(2, 1, 2, 'Parent', fig);
    plotTraceAxes(ax2, relGrid, yN, colNeg, 1.5, ...
        sprintf('Mouse: %s | Fos- Peri-Blink Mean Ca (n=%d cells)', id, sum(string(mouseNeg) == id & okNeg)), 9, 'normal');
    xlabel(ax2, 'Time from Blink Onset (sec)');
    appendImage(fig, savePath);
end

clf(fig);
ax1 = subplot(2, 1, 1, 'Parent', fig);
plotTraceAxes(ax1, relGrid, mean(tracePos, 1, 'omitnan')', colPos, 2.5, 'Overall Group Mean | Fos+ Peri-Blink Ca', 10, 'bold');
ax2 = subplot(2, 1, 2, 'Parent', fig);
plotTraceAxes(ax2, relGrid, mean(traceNeg, 1, 'omitnan')', colNeg, 2.5, 'Overall Group Mean | Fos- Peri-Blink Ca', 10, 'bold');
xlabel(ax2, 'Time from Blink Onset (sec)');
appendImage(fig, savePath);
end

function appendTimeCoursePage(binCenters, meanPos, semPos, meanNeg, semNeg, pHolm, pGroup, pInter, savePath, binWidth, aggLabel)
colPos = [0.85 0.1 0.1]; colNeg = [0.05 0.05 0.55];
fig = figure('Units', 'pixels', 'Position', [100 100 700 500], 'Color', 'w', 'Visible', 'off');
cleanup = onCleanup(@() close(fig));
ax = axes('Parent', fig, 'FontName', 'Arial', 'FontSize', 9, 'LineWidth', 0.75); hold(ax, 'on');

yAll = [meanPos + semPos; meanNeg + semNeg; meanPos - semPos; meanNeg - semNeg];
yAll = [yAll(~isnan(yAll)); 0];
yRange = max(yAll) - min(yAll); if yRange == 0, yRange = 1; end
yMin = min(yAll) - 0.05 * yRange;
yMax = max(yAll) + 0.20 * yRange;
xMax = max(binCenters);
patch(ax, [0 xMax xMax 0], [yMin yMin yMax yMax], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
text(ax, xMax/2, yMax - 0.05*yRange, 'Post-Blink', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'center');

errorbar(ax, binCenters, meanPos, semPos, '-s', 'Color', colPos, 'MarkerFaceColor', colPos, ...
    'MarkerEdgeColor', colPos, 'MarkerSize', 6, 'LineWidth', 1.3, 'CapSize', 3);
errorbar(ax, binCenters, meanNeg, semNeg, '-o', 'Color', colNeg, 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', colNeg, 'MarkerSize', 6, 'LineWidth', 1.3, 'CapSize', 3);
line(ax, [0 0], [yMin yMax], 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1);

for b = 1:numel(binCenters)
    p = pHolm(b);
    if isnan(p) || p >= 0.05, continue; end
    if p < 0.001, sig = '***'; elseif p < 0.01, sig = '**'; else, sig = '*'; end
    top = max([meanPos(b) + semPos(b), meanNeg(b) + semNeg(b)]);
    if isnan(top), continue; end
    text(ax, binCenters(b), top + 0.04 * yRange, sig, 'FontName', 'Arial', 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

xlabel(ax, 'Time Relative to Blink Onset (sec)', 'FontName', 'Arial', 'FontSize', 10);
ylabel(ax, 'Ca^{2+} Activity (DFF)', 'FontName', 'Arial', 'FontSize', 10);
xlim(ax, [min(binCenters) - 0.5, max(binCenters) + 0.5]);
xticks(ax, binCenters); ylim(ax, [yMin yMax]);

if ~isnan(pGroup) && ~isnan(pInter)
    sub = sprintf('Two-way LMM: Group p=%.3g, Group\\times{}Time p=%.3g (post-hoc per %g-sec bin, Holm-corrected)', pGroup, pInter, binWidth);
else
    sub = 'Two-way LMM: insufficient data';
end
title(ax, sprintf('Peri-Blink Ca^{2+} Activity Time Course (%g-sec bins, %s)', binWidth, aggLabel), 'FontName', 'Arial', 'FontSize', 11);
subtitle(ax, sub, 'FontName', 'Arial', 'FontSize', 8, 'Interpreter', 'tex');

hP = plot(ax, nan, nan, '-s', 'Color', colPos, 'MarkerFaceColor', colPos, 'MarkerEdgeColor', colPos);
hN = plot(ax, nan, nan, '-o', 'Color', colNeg, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colNeg);
legend([hP hN], {'Fos+', 'Fos-'}, 'Location', 'northeastoutside', 'Box', 'off', 'FontSize', 9);
set(ax, 'Box', 'off', 'TickDir', 'out');

appendImage(fig, savePath);
end

function T = appendMetaRow(T, analysisName, metricName, value, details)
if nargin < 4, value = NaN; end
if nargin < 5, details = ''; end
T = [T; table(string(analysisName), string(metricName), double(value), string(details), ...
    'VariableNames', T.Properties.VariableNames)];
end

function T = appendParamsToMeta(T, params)
desc = containers.Map();
desc('caTimeUnit')             = 'Time unit of the Ca data';
desc('fps')                    = 'Frame rate of the DLC video (Hz)';
desc('upperY_col')             = 'DLC: column index of upper eyelid Y coordinate';
desc('lowerY_col')             = 'DLC: column index of lower eyelid Y coordinate';
desc('win_sec_dee')            = 'One-sided window width of the moving median for DEE (sec)';
desc('f0Window')               = 'Time range for F0 calculation [start end] (sec)';
desc('f0Percentile')           = 'Percentile used for F0 calculation';
desc('analysisWindow')         = 'Time range for blink detection and Peri-Blink analysis [start end] (sec)';
desc('blinkDEEThreshold')      = 'Blink detection threshold (DEE_Inverted)';
desc('blinkDurationThreshold') = 'Minimum blink duration included in analysis (sec)';
desc('periBlinkWindowSec')     = 'Time range before/after blink onset to extract (sec)';
desc('periBlinkZeroAtOnset')   = 'Whether to zero each event at onset (0 sec)';
desc('timeCourseBinWidthSec')  = 'Bin width for time course graph and LMM (sec)';
desc('smoothWinSec_dee')       = 'Moving average window for DEE_Inverted (sec); disabled if <= 0';
desc('smoothWinSec_caDFF')     = 'Moving average window for Ca DFF (sec); disabled if <= 0';

fields = fieldnames(params);
for i = 1:numel(fields)
    f = fields{i}; v = params.(f);
    if isKey(desc, f), d = desc(f); else, d = ''; end
    if (isnumeric(v) || islogical(v)) && isscalar(v)
        T = appendMetaRow(T, 'Parameters', f, double(v), d);
    elseif isnumeric(v) || islogical(v)
        T = appendMetaRow(T, 'Parameters', f, NaN, sprintf('%s [value: %s]', d, mat2str(v)));
    else
        T = appendMetaRow(T, 'Parameters', f, NaN, sprintf('%s [value: %s]', d, char(string(v))));
    end
end
end

function T = appendMethodToMeta(T)
T = appendMetaRow(T, 'Analysis_Method', 'DFF', NaN, ...
    'For each cell, F0 is the f0Percentile-th percentile of Ca within f0Window, and DFF=(F-F0)/F0 is calculated within analysisWindow (optionally smoothed).');
T = appendMetaRow(T, 'Analysis_Method', 'Blink_Detection', NaN, ...
    'Eye opening is computed from the upper/lower eyelid Y coordinates in DLC. DEE=(x-moving median)/moving median is sign-inverted to give DEE_Inverted. Contiguous periods where DEE_Inverted >= blinkDEEThreshold are defined as blinks, and those with duration >= blinkDurationThreshold are analyzed.');
T = appendMetaRow(T, 'Analysis_Method', 'PeriBlink_Trace', NaN, ...
    'DFF within +/- periBlinkWindowSec sec around each blink onset (0 sec) is extracted and averaged per cell (all blinks of the same mouse). If periBlinkZeroAtOnset=true, the onset value of each event is subtracted before averaging. Events extending beyond the recording range are excluded.');
T = appendMetaRow(T, 'Analysis_Method', 'PeriBlink_TimeCourse', NaN, ...
    'DFF is averaged within bins of timeCourseBinWidthSec (PeriBlink_TimeCourse_Data sheet, 1 row = 1 cell). After an overall Two-way LMM (Value ~ Group*TimeBin + (1|Mouse) + (1|Mouse:Cell)), post-hoc LMM (Value ~ Group + (1|Mouse)) is run in each bin and corrected by Holm-Bonferroni. The last two TIFF pages show All Cells (cell-level mean +/- SEM) and Mouse-level (mouse-level mean +/- SEM); the statistical results are shared (cell-level LMM).');
end
