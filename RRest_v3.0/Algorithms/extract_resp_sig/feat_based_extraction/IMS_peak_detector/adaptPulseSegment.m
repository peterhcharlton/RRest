function [peaks,onsets,clipp] = adaptPulseSegment(y,Fs,annot)
%ADAPTPULSESEGMENT perform adaptive pulse segmentation and artifact detection 
%in ppg signals
%   [peaks,onsets,artif] = adaptPulseSegment(y,annot)
%
% Inputs:
%       y      vector, ppg signal [Lx1] or [1xL], in which L = length(signal)
%       Fs      scalar, sampling frequency in Hz
%       annot   vector, timestamps (in samples) with location of the peaks
%
% Outputs:
%       peaks   vector, locations of peaks (in samples)
%       onsets  vector, locations of onsets of the beats (in samples)
%       artif   vector, locations of peaks classified as artefacts
% 
% References:
%       Karlen et al., Adaptive Pulse Segmentation and Artifact Detection in 
%       Photoplethysmography for Mobile Applications, 34th Internat. Conf. IEEE-EMBS 2012
%       
% Written by Marco A. Pimentel

doOptimise = 1;
doPlot = 0;
if nargin < 3
    % no annotations are provided, therefore, no optimisation will take
    % place
    doOptimise = 0; 
end


% The algorithm in the paper is applied to signals sampled at 125 Hz...
% We do not resample our signal
%ys = resample(y,125,Fs);
%Fs = 125;
% if Fs ~= 300
%     ys = resample(y,300,Fs);
%     Fs = 300;
% else
    ys = y;
% end

% The paper is not clear about the selection of the range of search for m
% We define the range of m to be between [0.005 - 0.100] secs (5ms to 100ms)
% We define "m" in terms of samples
opt.bounds = 0.005:0.005:0.100;
opt.m = unique(ceil(opt.bounds*Fs));

opt.perf = zeros(length(opt.m),4); % store results of performance for each m

if doOptimise
    % Perform optimisation
    for i = 1 : length(opt.m)
        % Determine peaks and beat onsets
        [linez,linezSig] = pulseSegment(ys,Fs,opt.m(i));
        % Calculate performance of the peak detection
        opt.perf(i,:) = evalPerf(annot,linez(:,2));
    end
    
else
    % Do not perform optimization; fix m
    opt.m = 10;
    [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,opt.m);
end

if doPlot
    colData = {'g','y','r'};
    figure; 
    h(1) = subplot(211);
    plot(ys); hold on;
    for i = 1 : size(linez,1)
        %if linezSig(i) > -1
        plot(linez(i,:),ys(linez(i,:)),'-x','Color',colData{linezSig(i)+2});
        %end
    end
    
    h(2) = subplot(212);
    plot(ys,'g'); hold on;
    for i = 1 : size(peaks,1)
        plot(peaks(i,:),ys(peaks(i,:)),'-xr');
    end
    if ~isempty(artifs)
    for i = 1 : size(artifs,1)
        plot(artifs(i,:),ys(artifs(i,:)),'--^b');
    end
    end
    if ~isempty(clipp)
    for i = 1 : size(clipp,1)
        plot(clipp(i,:),ys(clipp(i,:)),'-om');
    end
    end
    linkaxes(h,'x');
    
end

% Correct for the downsmapling performed during the peak detection
onsets = peaks(:,1);
peaks  = peaks(:,2);
for i = 1 : size(peaks,1)
    [~,ind]  = min(ys(max([1 onsets(i)-opt.m]):min([length(ys) onsets(i)+opt.m])));
    onsets(i) = max([1 onsets(i)-opt.m]) + ind(1) - 1;
    [~,ind]  = max(ys(max([1 peaks(i)-opt.m]):min([length(ys) peaks(i)+opt.m])));
    peaks(i) = max([1 peaks(i)-opt.m]) + median(ind) - 1;
end

% Correct minimum value of onset of the beat
for i = 2 : length(onsets)
    [~,ind]   = min(ys(peaks(i-1):peaks(i)));
    onsets(i) = peaks(i-1) + ind - 1;
end

end

function [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,m)
% Perform pulse segmentation in ys given m
% Inputs:
%       ys      vector, with ppg signal
%       m       scalar, with the length of each line (read paper for details)
% 
% Outputs:
%       linez      2-column vector, with location of beat onsets and peaks
%       linezSig   vector, with label of the slope of each line segment
%                  1 - positive slope; -1 - negative slope; 0 - constant
% 

% split signal in different segments
nseg = floor(length(ys)/m);    % number of segments
% nseg = round(length(ys)/m);    % number of segments
% intialize loop variables
seg = 1;    % segment counter
z = 1;      % line counter 
segInLine = 1;  % line controler
linez = zeros(nseg,2); linez(1,:) = [1,m];
% slope of segment/line
a = zeros(nseg,1); a(1) = slope(ys,linez(1,:));
% "classify" line segment according to the slope
linezSig = zeros(nseg,1); linezSig(1) = sign(a(1));
% Start loop over segments
z = z + 1; seg = seg + 1;
for i = 2 : nseg    % loop over segments
    linez(z,:) = [(seg-1)*m+1 seg*m];
    try
        a(z) = slope(ys,linez(z,:));
    catch
        a = 1;
    end
    linezSig(z) = sign(a(z));
    if sign(a(z)) == sign(a(z-1))
        linez(z-1,:) = [(seg-1-segInLine)*m+1 seg*m];
        seg = seg + 1;
        segInLine = segInLine + 1;
    else
        z = z + 1;
        seg = seg + 1;
        segInLine = 1;
    end
end

% remove extra spaces created in output variables
linezSig(sum(linez,2)==0,:) = [];
linez(sum(linez,2)==0,:) = [];

% Apply adaptive threshold algorithm
% For this algorithm to work, we need to first find a valide line segment 
% in order to intialize the thresholds! In order to this, we define a flag
% to control the intialization in the main loop
FOUND_L1 = 0;

% The algorithm includes the definition of 4 adaptation parameters
% We define the following adaptation parameters
% a =     | a_fast_low    a_fast_high |
%         | a_slow_low    a_slow_high |
% 
a = [0.5 1.6; ...
     0.6 2.0];
 
% Define fixed thresholds described in the paper
ThT  = 0.03 * Fs;    % Duration of the line
ThIB = 0.24 * Fs;    % Interbeat invertal (240 ms) 

% Define parameters used in the main loop
alpha = zeros(size(linez,1),1);
for i = 1 : size(linez,1)
    alpha(i) = slope(ys,linez(i,:));   % slopes of line segments
end
theta = diff(ys(linez),[],2);
durat = diff(linez,[],2);       % duration of line segments (in samples)

% remove lines that do not have the necessary duration
linez(durat<ThT,:) = [];
theta(durat<ThT,:) = [];
alpha(durat<ThT,:) = [];
horiz = horizontalLine(ys,linez,Fs);

FLAG = 0;
artifs = []; clipp = [];
% Select window for detect firs peaks!
wind = theta(theta>0);
try 
    wind = wind(1:10);
catch
    wind = wind;
end
ThAlow  = prctile(wind,95)*0.6;
ThAhigh = prctile(wind,95)*1.8;
peaks = [];
for z = 1 : size(linez,1)-1   % loop over line segments
    if FOUND_L1
        if alpha(z) > 0 && ... % slope must be positive
                alpha(z-1) ~= 0 && ...  % peaks before or after clipping are artefactual
                alpha(z+1) ~= 0
            if theta(z) >= ThAlow && theta(z) <= ThAhigh && ...
                    linez(z,2) >= peaks(end,2) + ThIB
                ThAlow  = (ThAlow + theta(z)*a(2,1))/2;
                ThAhigh = theta(z) * a(2,2);
                FLAG = 0;
                currTheta = [currTheta; theta(z)];
                peaks = [peaks; linez(z,:)];
            else
                if FLAG > 0
                    ThAlow  = (ThAlow + min(currTheta(max([1 end-4]):end))*a(1,1))/2;
                    ThAhigh = max(currTheta(max([1 end-4]):end)) * a(1,2);
                    %ThAlow  = (ThAlow + theta(z)*a(1,1))/2;
                    %ThAhigh = theta(z) * a(1,2);
                end
                FLAG = FLAG + 1;
                artifs = [artifs; linez(z,:)];
            end
        elseif theta(z) > 0 && ... 
                ((theta(z-1) ~= 0 || horiz(z-1) ~= 0) && ...
                (theta(z+1) ~= 0 || horiz(z+1) ~= 0))
            if theta(z) >= ThAlow && theta(z) <= ThAhigh && ...
                    linez(z,2) >= peaks(end,2) + ThIB
                ThAlow  = (ThAlow + theta(z)*a(2,1))/2;
                ThAhigh = theta(z) * a(2,2);
                FLAG = 0;
                currTheta = [currTheta; theta(z)];
                peaks = [peaks; linez(z,:)];
            else
                if FLAG > 0
                    %ThAlow  = (ThAlow + currTheta*a(1,1))/2;
                    %ThAhigh = currTheta * a(1,2);
                    ThAlow  = (ThAlow + min(currTheta(max([1 end-4]):end))*a(1,1))/2;
                    ThAhigh = max(currTheta(max([1 end-4]):end)) * a(1,2);
                    %ThAlow  = (ThAlow + theta(z)*a(1,1))/2;
                    %ThAhigh = theta(z) * a(1,2);
                end
                FLAG = FLAG + 1;
                artifs = [artifs; linez(z,:)];
            end
        elseif theta(z) == 0 && horiz(z) == 0
            artifs  = [artifs; linez(z,:)];
            clipp   = [clipp; linez(z,:)];
        end 
    else
        if alpha(z) > 0 && durat(z) >= ThT && ...
                theta(z) >= ThAlow && theta(z) <= ThAhigh 
            FOUND_L1 = 1;
            ThAlow  = theta(z)*0.5;
            ThAhigh = theta(z)*2.0;
            peaks = linez(z,:);    % loaction of onsets and peaks
            currTheta = theta(z);
        end
    end
end

end

function out = horizontalLine(ys,linez,Fs)
% Get horizontal lines from signal given linez
out = zeros(size(linez,1),1);
for i = 1 : size(linez,1)
    out(i) = median(abs(diff(ys(linez(i,1):linez(i,2)))));
    % check duration of the peaks
    if out(i) == 0 && diff(linez(i,:)) <= 0.200*Fs
        out(i) = 0.1;
    end
end

end

function out = slope(ys,interv)
start = interv(1); stop = interv(2);
out = sum(diff(ys([start:stop])))/(stop-start);
%out = median(gradient(ys(start:stop)));
end