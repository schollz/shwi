function [peaks] = getPeaksWFFT(x,y,minPeakDistance,minPeakHeight,bestLength)

% x=data(:,1);
% y=data(:,2);
% Define vectors
forceData=y';
extensionData=x;
x=size(extensionData);
X=forceData;

%% Determine the window length, i.e. find datapoints/nm *roughly*
conversion=sum(extensionData>10)/(max(extensionData)-10);
windowLength = round(bestLength*conversion);
stepWindow=  round(conversion*0.5);
% Get length of data vector
dataLength = length(X);
windowedDataLength = ceil(dataLength/stepWindow);
if (windowLength/2 == round(windowLength/2))
	fftLength = (windowLength/2) + 1;
else
	fftLength = ceil(windowLength/2);
end

% Initialize vector for holding windowed FT
wfftSignal = zeros(windowedDataLength,fftLength);
X = [zeros(1,windowLength/2) X zeros(1,windowLength/2)]; % padding the signal with zeros
ii = 0;
% Loop through each window, taking FT
windowedSeries = 1:stepWindow:dataLength;
for iii=1:length(windowedSeries)
    i=windowedSeries(iii);
	ii = ii + 1;
	dataFrame = X(i:(i + windowLength - 1));
	fftValues = abs(fft(dataFrame, windowLength));
	wfftSignal(ii,:) = fftValues(1:fftLength);
end


% Sum up the odd coefficients from the FT
coefficientSum=zeros(length(wfftSignal(:,1)),1);
for i=1:2:size(wfftSignal,2)
    coefficientSum=coefficientSum+wfftSignal(:,i);
end
coefficientSum = coefficientSum/windowLength;

% Two cases: Signal ends in rupture or cyclic
% For rupture, simply baseline FT sum by taking average of final points
% (optimal)
% For cyclic, baseline FT by taking average of minimum of points
foo = sort(coefficientSum);
if (std(coefficientSum(end-50:end)) < 20)
    coefficientSum = coefficientSum-mean(coefficientSum(end-50:end));
else
    coefficientSum = coefficientSum-mean(foo(1:30));
end

PeakSig = (coefficientSum(1:end));

% Find peaks using the protein database parameters
[pks,locs] = findpeaks(PeakSig,'SORTSTR','descend','MinPeakProminence',5,'minpeakdistance',minPeakDistance,'minpeakheight',minPeakHeight);



expandingWindow = 150;
xPeaks = [];
yPeaks = [];
lastX = 0;
for j=1:length(locs)
    posVal = windowedSeries(locs(j));
    if posVal>expandingWindow && posVal<5000-expandingWindow
        areaForces = forceData(posVal-expandingWindow:posVal+expandingWindow);
        [maxval,maxind] = max(areaForces);
        [minval,minind] = min(areaForces);
        posVal = posVal - expandingWindow + maxind-1;
        ext = extensionData(posVal);
        xPeaks = [xPeaks; extensionData(posVal)];
        yPeaks = [yPeaks; forceData(posVal)];
    end
end


peaks = sortrows([unique(xPeaks),unique(yPeaks)]);
if length(peaks)>0
    peaks2 = [];
    clear diff;
    lcs = diff([peaks(1,1); peaks(:,1)]');
    for i=1:length(lcs)
        if lcs(i)<150
            peaks2 = [peaks2; peaks(i,:)];
        end
    end
    peaks = peaks2;
%
%     t = flipud(peaks);
%     for i=1:size(t,1)-1
%         if t(i,2)>50
%             break
%         end
%     end
%     t = t(i:end,:);
%     peaks = flipud(t);

end

%plot(extensionData,forceData,'-',peaks(:,1),peaks(:,2),'kv','markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0])

end
