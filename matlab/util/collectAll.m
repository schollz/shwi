
filelist = getAllFiles('./');
for i=1:length(filelist)
   if  length(findstr('.afm',char(filelist(i))))>0
       data=importdata(char(filelist(i)));
       [xPeaks,yPeaks] = analyzeCurve(data(:,1),data(:,2),10,20,5);
       plot(data(:,1),data(:,2),'-',xPeaks,yPeaks,'kv','markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0])
       pause(1)
       close all;
   end
end