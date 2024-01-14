%%%%read viral load data
tumordata=readtable('UT_NSC_Lung_Cancer');

TT(:,5) = tumordata.P5T;
TV(:,5)= tumordata.P5V;
%convert date to mm^3
TV(:,5)=TV(:,5)*10^3;
TT(:,4) = tumordata.P4T;
TV(:,4)= tumordata.P4V;
TV(:,4)=TV(:,4)*10^3;
TT(:,3) = tumordata.P3T;
TV(:,3)= tumordata.P3V;
TV(:,3)=TV(:,3)*10^3;
TT(:,2) = tumordata.P2T;
TV(:,2)= tumordata.P2V;
TV(:,2)=TV(:,2)*10^3;
TT(:,1) = tumordata.P1T;
TV(:,1)= tumordata.P1V;
TV(:,1)=TV(:,1)*10^3;