
% v 01 04.04.2023
% read methylation files 
%
%

%meta = readtable('Meth samples - human, apes, mammoth.xlsx','Sheet');
meta= readtable('Meth samples - human, apes, mammoth.xlsx','Sheet','Ancient');

% delete empty rows 

empty_rows = cellfun('isempty', meta.Sample);
meta=meta(~empty_rows ,:);

%loc=cellfun('isempty', meta.Sample );

%meta=  rmmissing(meta);

% 'Range','A1:E11'
% a = 1;
% while(a<24)
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1116.mat');
%     samp0 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5319.mat');
%     samp1 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4532.mat');
%     samp2 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2862.mat');
%     samp3 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4529.mat');
%     samp4 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_3758.mat');
%     samp5 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5725.mat');
%     samp6 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2861.mat');
%     samp7 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_3957.mat');
%     samp8 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_3565.mat');
%     samp9 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4792.mat');
%     samp10 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4775.mat');
%     samp11 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4315.mat');
%     samp12 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1053.mat');
%     samp13 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_7421.mat');
%     samp14 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_3255.mat');
%     samp15 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5835.mat');
%     samp16 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2514.mat');
%     samp17 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5748.mat');
%     samp18 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1633.mat');
%     samp19 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5950.mat');
%     samp20 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5838.mat');
%     samp21 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5742.mat');
%     samp22 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5743.mat');
%     samp23 = simu_fas.getmethylation(a);
%  %   load('C:\Users\amitai\Google Drive\p_CpGs\simu_2105.mat'); % no simu
%  %   pasturalist
%   %  samp24 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2520.mat');
%     samp25 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2935.mat');
%     samp26 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2978.mat');
%     samp27 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2980.mat');
%     samp28 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_3133.mat');
%     samp29 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4634.mat');
%     samp30 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1965.mat');
%     samp31 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1631.mat');
%     samp32 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1632.mat');
%     samp33 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1961.mat');
%     samp34 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_LBK.mat'); 
%     samp35 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1496.mat');
%     samp36 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5077.mat');
%     samp37 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1962.mat');
%     samp38 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4438.mat');
%     samp39 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_I15.mat');
%     samp40 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2134.mat');
%     samp41 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1507.mat');
%     samp42 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4878.mat');
%     samp43 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4432.mat');
%     samp44 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4873.mat');
%     samp45 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_Los.mat');
%     samp46 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_LaB.mat'); 
%     samp47 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4596.mat');
%     samp48 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5233.mat');
%     samp49 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_0708.mat');
%     samp50 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1960.mat');
%     samp51 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4914.mat');
%     samp52 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4875.mat');
%     samp53 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4877.mat');
%     samp54 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2139.mat');
%     samp55 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_SF1.mat'); 
%     samp56 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1734.mat');
%     samp57 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5236.mat');
%     samp58 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5235.mat');
%     samp59 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_M45.mat'); 
%     samp60 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_Ust.mat'); 
%     samp61 = simu_fas.getmethylation(a);
% 
%     fn = strcat('chr_',num2str(a),'_simu_united.csv') ;
% 
%     T = [samp0 samp1 samp2 samp3 samp4 samp5 samp6 samp7 samp8 samp9 samp10 samp11 samp12 samp13 samp14 samp15 samp16 samp17 samp18 samp19 samp20 samp21 samp22 samp23 samp25 samp26 samp27 samp28 samp29 samp30 samp31 samp32 samp33 samp34 samp35 samp36 samp37 samp38 samp39 samp40 samp41 samp42 samp43 samp44 samp45 samp46 samp47 samp48 samp49 samp50 samp51 samp52 samp53 samp54 samp55 samp56 samp57 samp58 samp59 samp60 samp61 ];
% 
%     csvwrite(fn, T);
% 
%     a = a+1;
% end


% a = 1;
% while(a<24)
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1116.mat');
%     samp0 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5319.mat');
%     samp1 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4532.mat');
%     samp2 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2862.mat');
%     samp3 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4529.mat');
%     samp4 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_3758.mat');
%     samp5 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5725.mat');
%     samp6 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2861.mat');
%     samp7 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_3957.mat');
%     samp8 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_3565.mat');
%     samp9 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4792.mat');
%     samp10 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4775.mat');
%     samp11 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4315.mat');
%     samp12 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1053.mat');
%     samp13 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_7421.mat');
%     samp14 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_3255.mat');
%     samp15 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5835.mat');
%     samp16 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2514.mat');
%     samp17 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5748.mat');
%     samp18 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1633.mat');
%     samp19 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5950.mat');
%     samp20 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5838.mat');
%     samp21 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5742.mat');
%     samp22 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5743.mat');
%     samp23 = simu_fas.getmethylation(a);
%  %   load('C:\Users\amitai\Google Drive\p_CpGs\simu_2105.mat'); % no simu
%  %   pasturalist
%   %  samp24 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2520.mat');
%     samp25 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2935.mat');
%     samp26 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2978.mat');
%     samp27 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2980.mat');
%     samp28 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_3133.mat');
%     samp29 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4634.mat');
%     samp30 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1965.mat');
%     samp31 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1631.mat');
%     samp32 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1632.mat');
%     samp33 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1961.mat');
%     samp34 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_LBK.mat'); 
%     samp35 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1496.mat');
%     samp36 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5077.mat');
%     samp37 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1962.mat');
%     samp38 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4438.mat');
%     samp39 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_I15.mat');
%     samp40 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2134.mat');
%     samp41 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1507.mat');
%     samp42 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4878.mat');
%     samp43 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4432.mat');
%     samp44 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4873.mat');
%     samp45 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_Los.mat');
%     samp46 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_LaB.mat'); 
%     samp47 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4596.mat');
%     samp48 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5233.mat');
%     samp49 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_0708.mat');
%     samp50 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1960.mat');
%     samp51 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4914.mat');
%     samp52 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4875.mat');
%     samp53 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_4877.mat');
%     samp54 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_2139.mat');
%     samp55 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_SF1.mat'); 
%     samp56 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_1734.mat');
%     samp57 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5236.mat');
%     samp58 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_5235.mat');
%     samp59 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_M45.mat'); 
%     samp60 = simu_fas.getmethylation(a);
%     load('C:\Users\amitai\Google Drive\p_CpGs\simu_Ust.mat'); 
%     samp61 = simu_fas.getmethylation(a);
% 
%     fn = strcat('chr_',num2str(a),'_simu_united.csv') ;
% 
%     T = [samp0 samp1 samp2 samp3 samp4 samp5 samp6 samp7 samp8 samp9 samp10 samp11 samp12 samp13 samp14 samp15 samp16 samp17 samp18 samp19 samp20 samp21 samp22 samp23 samp25 samp26 samp27 samp28 samp29 samp30 samp31 samp32 samp33 samp34 samp35 samp36 samp37 samp38 samp39 samp40 samp41 samp42 samp43 samp44 samp45 samp46 samp47 samp48 samp49 samp50 samp51 samp52 samp53 samp54 samp55 samp56 samp57 samp58 samp59 samp60 samp61 ];
% 
%     csvwrite(fn, T);
% 
%     a = a+1;
% end
% 
% 
