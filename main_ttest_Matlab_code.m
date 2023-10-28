clc
clear all
% Sample data (replace these arrays with your actual data)
load Features_before.mat;
load Features_during.mat;
data_before=Features_before(:,1,19)';
data_after=Features_during(:,1,19)';
k=data_after-data_before;
ttest_results=zeros(36,19);
for i=1:36
    for j=1:19

data_before=Features_before(:,j,i)';
data_after=Features_during(:,j,i)';
% Perform paired t-test
[h, p, ci, stats] = ttest(data_after -data_before);


ttest_results(i,j)=p;
   
 end
end

