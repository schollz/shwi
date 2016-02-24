bins = 0:4:500
[n1] = hist(lcs,bins)
[n2] = hist(lcs2,bins)
plot(bins,n1,bins,n2)

n2e = [zeros(1,100) n2 zeros(1,100)];
minResidual = 10000000000;
bestI = 0;
for i=0:200
    newn2 = n2e((1:length(n2))+i);
    subplot(2,1,1)
   plot(bins,n1,bins,newn2)
   subplot(2,1,2)
   metric = sum((newn2-n1).^2);
   plot(bins,(newn2-n1).^2)
   title(sprintf('sum sq: %2.1f',metric))
   pause(0.005)
   if metric < minResidual
       minResidual = metric;
       bestI = i;
   end
end
minResidual
bestI

subplot(2,1,1)
plot(data(:,1),data(:,2),data2(:,1),data2(:,2))
