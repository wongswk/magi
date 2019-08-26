function [meanFactor, sdFactor] = getFrequencyBasedPrior(x)

  z = fft(x);
  zmod = abs(z);
  zmod(1) = [];
  zmodEffective = zmod(1:(size(zmod,1)/2));
%   if(showplot) {
%     names(zmodEffective) <- 1:length(zmodEffective)
%     barplot(zmodEffective)
%     boxplot(zmodEffective)
%   }
  outliers = isoutlier(zmodEffective);
  if (sum(outliers)==0)
    [foo, freq] = max(zmodEffective);
  else
    outliers = isoutlier(zmodEffective) & zmodEffective > median(zmodEffective);
    %outliers = outliers[outliers > median(zmodEffective)]
    %whichOutliers <- which(zmodEffective %in% outliers)
    freq = find(outliers, 1, 'last');
  end
  meanFactor = 0.5 / freq;
  sdFactor = (1 - meanFactor) / 3;
  %c(meanFactor=meanFactor, sdFactor=sdFactor)
end