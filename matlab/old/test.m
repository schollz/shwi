%% ANALYSIS PARAMETERS
% Parameters for peak detection
minPeakDistance = 20; % Approximately min distance between peaks you'd expect
minPeakHeight = 50; % Approximately the min rupture force of peak
bestLength = 5; % Window size (nm). 5nm is recommended
% Parameters for WLC fit
P = 0.4; % Persitence length
kT = 4.1; % pN / nm. Don't change unless you are doing something with temperature


%% FILE COLLECTION
% Files can be kept in folders relating to the date/conditions
% Individual files are simply two columns of the extension and force data
% If you are recording displacement and force data, you will need to
% transform your data
fileList =  getAllFiles('8i27');

%% DATA ANALYSIS
numData = 0;
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:length(fileList)
    a = importdata(char(fileList(i)));
    x = a(:,1);
    y = a(:,2);
    data.originalX = x;
    data.originalY = y;
    
    % Get baseline tries to format the curve so that it starts at (0,0)
    [baselines] = getBaseline(x,y);
    x = x - baselines(1);
    y = y - baselines(2);
    data.x = x;
    data.y = y;
    
    % Peak detection
    [peaks] = analyzeCurve(x,y,minPeakDistance,minPeakHeight,bestLength);
    data.peaks = peaks;
    
    % Get contour-length increments
    Ls = [];
    for j=1:size(peaks,1)
        F = peaks(j,2);
        ext = peaks(j,1);
        myfun = @(L) F*P/kT-1/(4*(1-ext/L).^2)+1/4-ext/L;
        L0=fzero(myfun,round(ext)*1.5); 
        Ls = [Ls; L0];
    end
    dLs = diff(Ls);
    data.Ls = [Ls; 0];

    % Save data
    numData = numData + 1;
    datas{numData} = data;
 
    % Plot the examples (comment this out if you don't want it)
    plot(x,y,'-',peaks(:,1),peaks(:,2),'o','MarkerSize',20) 
    hold on;
    for j=1:length(dLs)
        text(mean([peaks(j,1) peaks(j+1,1)]),mean([peaks(j,2) peaks(j+1,2)]),sprintf('%2.1f',dLs(j)));
    end
    for j=1:length(Ls)
         plot(0:round(Ls(j)),WLCcurve(P,Ls(j),0:round(Ls(j))),'r');
    end
    hold off;
    axis([min(x) max(x) min(y) 1.5*max(y)])
    pause(1)
    
end



%% PLOTTING
% Outliers are probably due to problems with automatic peak detection
Fs = [];
Ls = [];
for i=1:length(datas)
    Fs = [Fs; datas{i}.peaks(1:end-1,2)];
    Ls = [Ls; diff(datas{i}.Ls(1:end-1))];
end
plot(Ls,Fs,'o')
xlabel('Contour-length increment [nm]')
ylabel('Force [pN]')

