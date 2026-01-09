close all; clear; clc;

load('G:\.shortcut-targets-by-id\1m9bxSOw_XjRwQOO_UZfscrFFu-hDWggr\ME310 Lab Friday Group 3\Project\24Nov\CalibrationVDAmp\calibrationfit.mat');

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = [24, 10023];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = "x";
opts.VariableTypes = "double";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

x = 0:9999;
data0 = table2array(readtable('VDAmpData0.csv', opts));
data10 = table2array(readtable('VDAmpData10.csv', opts));
data15 = table2array(readtable('VDAmpData15.csv', opts));
data20 = table2array(readtable('VDAmpData20.csv', opts));
data25 = table2array(readtable('VDAmpData25.csv', opts));
data30 = table2array(readtable('VDAmpData30.csv', opts));
data35 = table2array(readtable('VDAmpData35.csv', opts));
data40 = table2array(readtable('VDAmpData40.csv', opts));
data45 = table2array(readtable('VDAmpData45.csv', opts));
data50 = table2array(readtable('VDAmpData50.csv', opts));
data55 = table2array(readtable('VDAmpData55.csv', opts));
data60 = table2array(readtable('VDAmpData60.csv', opts));
data65 = table2array(readtable('VDAmpData65.csv', opts));
data70 = table2array(readtable('VDAmpData70.csv', opts));
data75 = table2array(readtable('VDAmpData75.csv', opts));
data80 = table2array(readtable('VDAmpData80.csv', opts));
data85 = table2array(readtable('VDAmpData85.csv', opts));
data90 = table2array(readtable('VDAmpData90.csv', opts));
data95 = table2array(readtable('VDAmpData95.csv', opts));
data100 = table2array(readtable('VDAmpData100.csv', opts));

%Combine all data to one array
combinedData = [data0 data10 data15 data20 data25 data30 data35 data40 data45 data50 data55 data60 data65 data70 data75 data80 data85 data90 data95 data100];

filteredData = combinedData(:,2:end);

%% lol why
%for i = 1:19
    %figure
    %plot(filteredData(:,i))
    %title(sprintf('Filtered Data for Trial %d', i));
    %xlabel('Samples');
    %ylabel('Amplitude');
    %grid on;
%end

%% Aamir's shit
% OMIT OR EDIT TRIAL 8 -- AVG PERIOD IS 146 (Trial 7 is 254 and Trial 9 is 183)

numTrials = size(filteredData,2);
periods = zeros(1, numTrials);   % store period for each trial
avgAmplitudes = zeros(1, numTrials);  % store average amplitude for each trial

for col = 1:numTrials
    
    % --- Extract column ---
    x = filteredData(filteredData(:,col) > 3.5, col);
    x = (x./a_1) - (a_0/10);
    
    % --- Quick frequency estimate (before heavy filtering) ---
    N = length(x);
    Y = fft(x);
    P2 = abs(Y/N);
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = (0:(floor(N/2)))/(N);  % Normalized frequency
    
    [~, maxIdx] = max(P1(2:end));  % Skip DC component
    estFreq = f(maxIdx + 1);
    
    % --- Adaptive preprocessing based on frequency ---
    if estFreq > 0.15  % High frequency (>~3 Hz)
        % Skip median filter entirely, use light moving average
        windowSize = 3;
        x1 = movmean(x, windowSize);
    else  % Low-mid frequency
        % Median filter is fine for slower signals
        x1 = medfilt1(x, 3);
    end
    
    % --- Adaptive filtering based on estimated frequency ---
    % For high frequencies, use minimal or no low-pass filtering
    if estFreq > 0.15  % High frequency (>3 Hz at typical sample rates)
        % Skip Butterworth entirely - just use the lightly smoothed data
        filteredx = x1;
    else
        % Normal adaptive filtering for low-mid frequencies
        cutoff = min(4 * estFreq, 0.45);
        cutoff = max(cutoff, 0.05);
        [b,a] = butter(2, cutoff);
        filteredx = filtfilt(b, a, x1);
    end
    
    % --- Zero-crossing detection approach ---
    % Remove DC offset
    mid = mean(filteredx);
    centered = filteredx - mid;
    
    % Find all sign changes (zero crossings)
    signChanges = diff(sign(centered));
    
    % Rising crossings: -1 or -2 -> +1 or +2 (signChanges == 2 or signChanges > 0)
    % Falling crossings: +1 or +2 -> -1 or -2 (signChanges == -2 or signChanges < 0)
    risingCrossings = find(signChanges > 0);
    fallingCrossings = find(signChanges < 0);
    
    % Combine and sort all crossings
    allCrossings = sort([risingCrossings; fallingCrossings]);
    
    % --- Build square wave ---
    sq = zeros(size(filteredx));
    
    % Determine initial state based on first sample
    if centered(1) >= 0
        currentState = 1;
    else
        currentState = 0;
    end
    
    sq(1:allCrossings(1)) = currentState;
    
    % Fill in square wave between crossings
    for i = 1:length(allCrossings)-1
        currentState = 1 - currentState;  % Toggle
        sq(allCrossings(i)+1:allCrossings(i+1)) = currentState;
    end
    
    % Fill to end
    currentState = 1 - currentState;
    sq(allCrossings(end)+1:end) = currentState;
    
    % --- Use FFT on square wave to find dominant frequency ---
    Nsq = length(sq);
    Ysq = fft(sq);
    P2sq = abs(Ysq/Nsq);
    P1sq = P2sq(1:floor(Nsq/2)+1);
    P1sq(2:end-1) = 2*P1sq(2:end-1);
    fsq = (0:(floor(Nsq/2)))/(Nsq);  % Normalized frequency
    
    % Find dominant frequency (skip DC component)
    [~, maxIdxSq] = max(P1sq(2:end));
    dominantFreq = fsq(maxIdxSq + 1);
    dominantPeriod = 1 / dominantFreq;
    
    % --- Reconstruct clean square wave at dominant frequency ---
    % Find phase by looking at first rising edge in noisy square wave
    sqDiff = diff([0; sq(:)]);
    noisyRiseIdx = find(sqDiff > 0);
    
    if ~isempty(noisyRiseIdx)
        % Use first rising edge as phase reference
        phaseOffset = noisyRiseIdx(1);
    else
        phaseOffset = 0;
    end
    
    % Generate clean square wave
    cleanSq = zeros(size(sq));
    t = (0:length(cleanSq)-1)';
    cleanSq = square(2*pi*dominantFreq*(t - phaseOffset)) > 0;
    cleanSq = double(cleanSq);
    
    % --- Identify rising edges in CLEAN square wave ---
    cleanSqDiff = diff([0; cleanSq(:)]);
    riseIdx = find(cleanSqDiff > 0);
    
    % --- Period calculation excluding first rising edge ---
    if length(riseIdx) >= 3  % Need at least 3 rising edges (skip first, need 2+ for intervals)
        % Skip first rising edge, calculate intervals from remaining
        riseIntervals = diff(riseIdx(2:end));
        periods(col) = mean(riseIntervals);
    elseif length(riseIdx) == 2
        % Only 2 rising edges - just use the one interval
        periods(col) = riseIdx(2) - riseIdx(1);
    else
        % Fallback to FFT-derived period
        periods(col) = dominantPeriod;
    end
    
    % --- Amplitude analysis: Iterative peak finding ---
    allPeaks = [];
    allPeakLocs = [];
    minPeakDistance = round(periods(col) * 0.4);  % Peaks should be ~half period apart minimum
    
    % Start with high prominence to find obvious peaks
    prominenceThreshold = std(filteredx) * 0.5;
    
    % Iteratively lower prominence to find more peaks
    maxIterations = 10;
    for iter = 1:maxIterations
        [pks, locs] = findpeaks(filteredx, 'MinPeakProminence', prominenceThreshold, ...
                                'MinPeakDistance', minPeakDistance);
        
        if ~isempty(pks)
            % Add new peaks
            allPeaks = [allPeaks; pks];
            allPeakLocs = [allPeakLocs; locs];
            
            % Check consistency: calculate coefficient of variation
            if length(allPeaks) >= 5
                peakStd = std(allPeaks);
                peakMean = mean(allPeaks);
                cv = peakStd / peakMean;  % Coefficient of variation
                
                % If peaks are consistent (CV < 0.15), we're done
                if cv < 0.15
                    break;
                end
            end
        end
        
        % Lower prominence threshold for next iteration
        prominenceThreshold = prominenceThreshold * 0.7;
        
        % Stop if threshold gets too low
        if prominenceThreshold < std(filteredx) * 0.05
            break;
        end
    end
    
    % Remove duplicate peaks (same location found in multiple iterations)
    [allPeakLocs, uniqueIdx] = unique(allPeakLocs);
    allPeaks = allPeaks(uniqueIdx);
    
    % Remove outlier peaks using MAD
    if length(allPeaks) >= 3
        medianPeak = median(allPeaks);
        mad = median(abs(allPeaks - medianPeak));
        
        if mad > 0
            % Keep peaks within 2.5 MAD of median (catches spikes)
            validMask = abs(allPeaks - medianPeak) < 2.5 * mad;
            validPeaks = allPeaks(validMask);
            validPeakLocs = allPeakLocs(validMask);
        else
            validPeaks = allPeaks;
            validPeakLocs = allPeakLocs;
        end
        
        % Calculate amplitude as distance from midline
        if ~isempty(validPeaks)
            peakAmplitudes = abs(validPeaks - mid);
            avgAmplitude = mean(peakAmplitudes);
            avgAmplitudes(col) = avgAmplitude;
        else
            avgAmplitude = NaN;
            avgAmplitudes(col) = NaN;
        end
    else
        validPeaks = allPeaks;
        validPeakLocs = allPeakLocs;
        if ~isempty(allPeaks)
            peakAmplitudes = abs(allPeaks - mid);
            avgAmplitude = mean(peakAmplitudes);
            avgAmplitudes(col) = avgAmplitude;
        else
            avgAmplitude = NaN;
            avgAmplitudes(col) = NaN;
        end
    end
    
    % --- Reconstruct thresholds for visualization ---
    stdVal = std(centered);
    highT = mid + 0.02 * stdVal;
    lowT = mid - 0.02 * stdVal;
    
    % --- Plotting ---
    figure(1); clf;
    
    subplot(2,1,1);
    plot(filteredx,'LineWidth',1.2,'Color',[0 0.4470 0.7410]); hold on;
    yline(mid,'--r','LineWidth',1.5);
    yline(highT,'--g','LineWidth',1);
    yline(lowT,'--g','LineWidth',1);
    % Plot all peaks found
    if ~isempty(allPeakLocs)
        plot(allPeakLocs, allPeaks, 'mx', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
    % Plot valid peaks (after outlier removal)
    if exist('validPeakLocs', 'var') && ~isempty(validPeakLocs)
        plot(validPeakLocs, validPeaks, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    end
    title(sprintf('Trial %d: Filtered Data (Est.Freq=%.4f)', col, estFreq));
    legend('Filtered Data','Mid','HighT','LowT','All Peaks','Valid Peaks','Location','best');
    grid on;
    
    subplot(2,1,2);
    plot(sq,'LineWidth',1.0,'Color',[0.7 0.7 0.7]); hold on;  % Noisy square wave in gray
    plot(cleanSq,'LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);  % Clean square wave
    plot(allCrossings, sq(allCrossings), 'ko', 'MarkerSize', 4);
    if ~isempty(riseIdx)
        plot(riseIdx(1), cleanSq(riseIdx(1)), 'rx', 'MarkerSize', 10, 'LineWidth', 2);  % First (excluded)
        if length(riseIdx) > 1
            plot(riseIdx(2:end), cleanSq(riseIdx(2:end)), 'g^', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        end
    end
    ylim([-0.2 1.2]);
    if ~isnan(avgAmplitude)
        title(sprintf('Square Wave | Period: %.1f samples | Avg Amplitude: %.3f (%d peaks)', ...
                      periods(col), avgAmplitude, length(validPeaks)));
    else
        title(sprintf('Square Wave | Period: %.1f samples', periods(col)));
    end
    legend('Noisy SqWave', 'Clean SqWave (FFT)', 'Zero Crossings', 'First Rise (excluded)', 'Used Rising Edges','Location','best');
    grid on;
    
    fprintf('\nTrial %d:\n', col);
    fprintf('  Est. Frequency: %.4f\n', estFreq);
    fprintf('  Dominant SqWave Frequency (FFT): %.4f (Period: %.1f samples)\n', dominantFreq, dominantPeriod);
    fprintf('  Total rising edges: %d (first excluded from calc)\n', length(riseIdx));
    fprintf('  Average Period: %.3f samples\n', periods(col));
    if ~isnan(avgAmplitude)
        fprintf('  Total peaks found: %d\n', length(allPeaks));
        fprintf('  Valid peaks (after outlier removal): %d\n', length(validPeaks));
        fprintf('  Average Amplitude from midline: %.3f\n', avgAmplitude);
        if length(validPeaks) > 1
            fprintf('  Peak std: %.3f (CV: %.2f%%)\n', std(validPeaks), 100*std(validPeaks)/mean(validPeaks));
        end
    end
    input('Press ENTER to continue to next trial...');
end

fprintf('\nAll periods computed.\n\n');
disp(array2table([periods; avgAmplitudes], ...
     'VariableNames', cellstr(string(1:numTrials)), ...
     'RowNames', {'Period_samples', 'Amplitude_from_midline'}));

%% Amplitude Analysis - Find maximum in each cycle

numTrials = size(filteredData,2);
avgAmplitudes = zeros(1, numTrials);

for col = 1:numTrials
    
    % --- Re-extract and filter data (same as before) ---
    x = filteredData(filteredData(:,col) > 3.5, col);
    
    N = length(x);
    Y = fft(x);
    P2 = abs(Y/N);
    P1 = P2(1:floor(N/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = (0:(floor(N/2)))/(N);
    
    [~, maxIdx] = max(P1(2:end));
    estFreq = f(maxIdx + 1);
    
    if estFreq > 0.15
        windowSize = 3;
        x1 = movmean(x, windowSize);
    else
        x1 = medfilt1(x, 3);
    end
    
    if estFreq > 0.15
        filteredx = x1;
    else
        cutoff = min(4 * estFreq, 0.45);
        cutoff = max(cutoff, 0.05);
        [b,a] = butter(2, cutoff);
        filteredx = filtfilt(b, a, x1);
    end
    
    % --- Generate clean square wave (same as before) ---
    mid = mean(filteredx);
    centered = filteredx - mid;
    
    signChanges = diff(sign(centered));
    risingCrossings = find(signChanges > 0);
    fallingCrossings = find(signChanges < 0);
    allCrossings = sort([risingCrossings; fallingCrossings]);
    
    sq = zeros(size(filteredx));
    if centered(1) >= 0
        currentState = 1;
    else
        currentState = 0;
    end
    sq(1:allCrossings(1)) = currentState;
    for i = 1:length(allCrossings)-1
        currentState = 1 - currentState;
        sq(allCrossings(i)+1:allCrossings(i+1)) = currentState;
    end
    currentState = 1 - currentState;
    sq(allCrossings(end)+1:end) = currentState;
    
    Nsq = length(sq);
    Ysq = fft(sq);
    P2sq = abs(Ysq/Nsq);
    P1sq = P2sq(1:floor(Nsq/2)+1);
    P1sq(2:end-1) = 2*P1sq(2:end-1);
    fsq = (0:(floor(Nsq/2)))/(Nsq);
    
    [~, maxIdxSq] = max(P1sq(2:end));
    dominantFreq = fsq(maxIdxSq + 1);
    
    sqDiff = diff([0; sq(:)]);
    noisyRiseIdx = find(sqDiff > 0);
    
    if ~isempty(noisyRiseIdx)
        phaseOffset = noisyRiseIdx(1);
    else
        phaseOffset = 0;
    end
    
    cleanSq = zeros(size(sq));
    t = (0:length(cleanSq)-1)';
    cleanSq = square(2*pi*dominantFreq*(t - phaseOffset)) > 0;
    cleanSq = double(cleanSq);
    
    cleanSqDiff = diff([0; cleanSq(:)]);
    riseIdx = find(cleanSqDiff > 0);
    
    % --- Find maximum in each cycle ---
    cycleMaxima = [];
    
    if length(riseIdx) >= 2
        % Skip first rising edge, analyze complete cycles
        for i = 2:length(riseIdx)
            % Define cycle from this rising edge to next (or end)
            if i < length(riseIdx)
                cycleStart = riseIdx(i);
                cycleEnd = riseIdx(i+1) - 1;
            else
                cycleStart = riseIdx(i);
                cycleEnd = length(filteredx);
            end
            
            % Extract cycle data
            cycleData = filteredx(cycleStart:cycleEnd);
            
            % Remove spikes: smooth with moving median to remove outliers
            smoothedCycle = movmedian(cycleData, 5);
            
            % Find maximum in smoothed cycle
            cycleMax = max(smoothedCycle);
            cycleMaxima(end+1) = cycleMax;
        end
        
        % Average amplitude across all cycles
        avgAmplitudes(col) = mean(cycleMaxima) - mid;
    else
        avgAmplitudes(col) = NaN;
    end
    
    % --- Visualization ---
    figure(2); clf;
    
    subplot(2,1,1);
    plot(filteredx, 'b', 'LineWidth', 1); hold on;
    plot(riseIdx(2:end), filteredx(riseIdx(2:end)), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    yline(mean(cycleMaxima), '--r', 'LineWidth', 2);
    title(sprintf('Trial %d: Cycle Analysis', col));
    legend('Filtered Data', 'Cycle Start Points', sprintf('Avg Amplitude: %.2f', avgAmplitudes(col)));
    xlabel('Samples');
    ylabel('Amplitude');
    grid on;
    
    subplot(2,1,2);
    if ~isempty(cycleMaxima)
        plot(cycleMaxima, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
        yline(avgAmplitudes(col), '--r', 'LineWidth', 2);
        title(sprintf('Maximum per Cycle (Mean: %.2f, Std: %.2f)', ...
                      avgAmplitudes(col), std(cycleMaxima)));
        xlabel('Cycle Number');
        ylabel('Maximum Amplitude');
        grid on;
    end
    
    fprintf('\nTrial %d Amplitude Analysis:\n', col);
    fprintf('  Number of complete cycles: %d\n', length(cycleMaxima));
    fprintf('  Average amplitude: %.3f\n', avgAmplitudes(col));
    if ~isempty(cycleMaxima)
        fprintf('  Amplitude std dev: %.3f\n', std(cycleMaxima));
    end
    
    %input('Press ENTER to continue to next trial...');
end

fprintf('\n=== Amplitude Summary ===\n');
disp(array2table([periods' avgAmplitudes']', 'VariableNames', ...
     cellstr(string(1:numTrials)), 'RowNames', {'Period_samples', 'Avg_Amplitude'}));

%% okay ... now to do the rest of ts

%For M(w), farny said output/input, both of that in ours is in meters or
%length

%change volt to length
m_values = (avgAmplitudes./a_1)./0.03175;
m_values = m_values ./ m_values(1);



periods(8) = 210; % :)
experimentalFreq = (1000 ./ periods)*2*pi;









%% Get omega_n via m_max

m_max = max(m_values);

syms z a

% Define the equation
eqn = -z^2 + z == a;

% Solve for z
sol_z = solve(eqn, z);

sol_sub = subs(sol_z, a, (1/(4*m_max^2))); % where a = 1/(4*m_max^2); derive from eqn 82 of class eqn doc

z = sqrt(double(sol_sub));
zeta = z(z < 0.707); % aamir hw8: when solving quadratic for zeta, take zeta<0.707

omega_r = experimentalFreq(find(m_values == m_max, 1)); % find index in expFreq that matches that of m_max

omega_n = omega_r/sqrt(1-(2*zeta^2));


fprintf('Zeta is %.4f.\n', zeta)
fprintf('Omega_n is %.4f.\n', omega_n)


%% Do sum plotting shits

% make arbitrary oomega vector
w = linspace(0,50);

M1 = 1 ./ sqrt((1 - (w./omega_n).^2).^2 + (2*zeta*(w./omega_n)).^2);
phi1 = -atan2(2*zeta*(w./omega_n), 1 - (w./omega_n).^2);

% Define the error value
errorValue = 0.0797;

figure(3);
plot(experimentalFreq(1:12), m_values(1:12), 'bo', 'LineWidth', 1);
ylim([0, 5]);
xlabel('angular frequency, rad/s');
ylabel('M(omega)');
title('M(omega) v.s. angular frequency');
plot(experimentalFreq(1:12), m_values(1:12), 'bo', 'LineWidth', 1);
errorbar(experimentalFreq(1:12), m_values(1:12), errorValue*ones(size(m_values(1:12))), 'b', 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1);
ylim([0, 5]);
xlabel('angular frequency, rad/s');
ylabel('M(omega)');
title('M(omega) v.s. angular frequency');
grid on;
hold on;
plot(w, M1, 'r')
legend('Experimentally-Determined M', 'Theoretical M')

figure(4);
plot(w, phi1, 'r')
xlabel('angular frequency, rad/s');
ylabel('phase shift, rad');
title('Theoretical Phase Shift v.s. angular frequency')
legend('Theoretical Phase Lag, rad')


%% calibrate observed freq with dial

% filteredData ranges from 10 to 100 in intervals of 5 (19 long)
% expFreq is also 19 long

% Create a table for dial percentage and experimental frequency
dialPercentage = (10:5:100)'; % Column vector from 10 to 100 in intervals of 5
experimentalFrequency = experimentalFreq(1:length(dialPercentage)); % Match length

% Compute third column (Hz instead of rad/s)
frequencyHz = experimentalFrequency / (2*pi);

% Create the table
freqTable = table(dialPercentage, experimentalFrequency', frequencyHz', ...
    'VariableNames', {'Dial Value, %', 'Experimental Freq, rad/s', 'Experimental Freq, Hz'});

% Display the table as a figure
figure(5); clf;

t = uitable('Data', freqTable{:,:}, ...
    'ColumnName', freqTable.Properties.VariableNames, ...
    'RowName', []);

set(t, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Extract the first 11 rows of the first two columns from freqTable
xData = freqTable{1:11, 'Dial Value, %'};
yData = freqTable{1:11, 'Experimental Freq, rad/s'};

% Perform polynomial fitting (linear fit)
p = polyfit(xData, yData, 1);

% Generate fitted values for plotting
fittedY = polyval(p, xData);

% Plot fit
figure(6); clf;
plot((freqTable{:, 'Dial Value, %'}), (freqTable{:, 'Experimental Freq, rad/s'}), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
plot(xData, fittedY, 'r-', 'LineWidth', 1.5);
xlabel('Dial Value, %');
ylabel('Experimental Frequency, rad/s');
title('Polynomial Fit of Experimental Frequency vs. Dial Value');
legend('Experimental Data', sprintf('Fit: y = %.4fx + %.4f', p(1), p(2)), 'Location', 'best');
grid on;


