load('Z:\Matlab_Data\TA_convperc_0708.mat');
load('Z:\Matlab_Data\TA_convperc_1053.mat');
load('Z:\Matlab_Data\TA_convperc_1116.mat');
load('Z:\Matlab_Data\TA_convperc_1496.mat');
load('Z:\Matlab_Data\TA_convperc_1507.mat');
load('Z:\Matlab_Data\TA_convperc_1631.mat');
load('Z:\Matlab_Data\TA_convperc_2520.mat');
samp1 = TA_convperc_0708;
samp2 = TA_convperc_1053;
samp3 = TA_convperc_1116;
samp4 = TA_convperc_1496;
samp5 = TA_convperc_1507;
samp6 = TA_convperc_1631;
samp7 = TA_convperc_2520;


for k=1:size(samp1,2)
    fn = strcat('chr_',num2str(k),'_united.csv') ;
    T1 = samp1{1,k};
    T2 = samp2{1,k};
    T3 = samp1{1,k};
    T4 = samp2{1,k};
    T5 = samp1{1,k};
    T6 = samp2{1,k};
    T7 = samp2{1,k};
    T = [T1 T2 T3 T4 T5 T6 T7];
%    labels = {'samp_0708'; 'samp_1053'; 'samp_1116'; 'samp_1496'; 'samp_1507'; 'samp_1631'; 'samp_2520'};
 %   data = [labels;T];
    csvwrite( fn, T );
end



