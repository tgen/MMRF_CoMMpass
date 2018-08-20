function roll_med=roll_med(sample,med_size)
% roll_med    -   Calculates average for a sliding window for different
% window lengths
%
%    Usage:
%
%         bpd_mins= roll_ave(sample,[3,5,9,15,25,51])
%
%    Input:
%         med_size  - Average window sizes to calculate  
%         sample    - The data to permutate
%
%    Output:
%
%         bpd_mins      - A matrix with med_size columns, each containing the med_size(i) sliding window calculations 
%  Obs:  NaN's are put for the first med_size(i) values.  
%   Written by: David Craig
%   Last Edit: October 24, 2005

[m,n]=size(sample);
[p]=length(med_size);
for i=1:p
    med_size(i);
    med_mat=sample;
    for g=2:med_size(i);
        med_mat(:,g)=circshift(med_mat(:,g-1),1);
    end
    med_mat=med_mat';
    med_mat=median(med_mat);
    roll_med(:,i)=med_mat';
    for h=1:med_size(i)
        roll_med(h,i)=NaN;
    end
end
