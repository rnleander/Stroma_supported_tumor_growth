%%read mouse data%%
%%%%read tumor growth data
tumordata=readtable('Untreated_JIMT_1_Mice');

TT(:,10) = tumordata.M10T;
TV(:,10)= tumordata.M10V;
TT(:,9) = tumordata.M9T;
TV(:,9)= tumordata.M9V;
TT(:,8) = tumordata.M8T;
TV(:,8)= tumordata.M8V;
TT(:,7) = tumordata.M7T;
TV(:,7)= tumordata.M7V;
TT(:,6) = tumordata.M6T;
TV(:,6)= tumordata.M6V;
TT(:,5) = tumordata.M5T;
TV(:,5)= tumordata.M5V;
TT(:,4) = tumordata.M4T;
TV(:,4)= tumordata.M4V;
TT(:,3) = tumordata.M3T;
TV(:,3)= tumordata.M3V;
TT(:,2) = tumordata.M2T;
TV(:,2)= tumordata.M2V;
TT(:,1) = tumordata.M1T;
TV(:,1)= tumordata.M1V;
