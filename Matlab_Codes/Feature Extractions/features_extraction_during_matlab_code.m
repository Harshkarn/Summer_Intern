load beta_during_all.mat
Features=zeros(4,19,size(During,3));
for i=1:size(During,3)
    for j=1:19
        Features(1,j,i)=mad(During(:,j,i));
        Features(2,j,i)=mad(During(:,j,i),1);
        Features(3,j,i)=std(During(:,j,i),1);
     
       histogram_values=histcounts(During(:,j,i),'Normalization','probability');
       histogram_values(histogram_values==0)=[];
       entropy_value=-sum(histogram_values.*log2(histogram_values));
       Features(4,j,i)=entropy_value;
    end
end
