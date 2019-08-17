function filled = insertNaN(mydata, level)
  if (level==0)
    filled = mydata;
    return
  end
  newdata = mydata;
  newdata = sortrows(newdata);
  dummydata = newdata(2:size(newdata,1),:);
  dummydata(:,:) = NaN;
  dummydata(:,1) = (newdata(2:size(newdata,1),1) + newdata(1:(size(newdata,1)-1),1))/2;
  newdata = vertcat(newdata, dummydata);
  newdata = sortrows(newdata);
  filled = insertNaN(newdata, level-1);
