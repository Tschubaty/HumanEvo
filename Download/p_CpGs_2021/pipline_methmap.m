%% General parameters
sample = 'Vin';
g_name = 'HG';
sex = 'F';
Library = 'single';
Organism_Kind = 'Homo Sapiens';

if strcmp(Organism_Kind,'Homo Sapiens')
    gc = load('gc_CpGs.mat');
    gc = gc.gc_CpGs;
    chr_names = gc.get('chr_names');
end

%% SAM to counts object
filename = [Organism_spec,'.sam'];
if strcmp(Organism_kind,'Homo Sapiens')
    %Human (hg19/hg38)
    if strcmp(Organism_spec,'Otzi') %hg38
        chr_lengths = [248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415];
        load HumanHg38
        GenomeAssembly = HumanHg38;
        reference = 'hg38';
    else %hg19
        chr_lengths = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566];
        load HumanHg19
        GenomeAssembly = HumanHg19;
        reference = 'hg19';
    end
elseif strcmp(Organism_kind,'Mammoth')
    chr_lengths = [214701375,225582790,205080690,166642186,143278944,163261726,83540442,128409435,94531929,103745917,68864707,83603808,103893473,70853943,88640724,45784276,69133140,72444150,68937429,79033276,62548349,54904812,53701227,27123610,58263816,76958368,50221267,120050768];
    load loxAfr4
    GenomeAssembly = loxAfr4;
    reference = 'loxAfr4';
else
    error('unrecognized organism')
end

if strcmp(Organism_kind,'Mammoth')
    fname = ['\\132.64.65.160\David.Gokhman\',Organism_kind,'\Methylome\Asian-African_fixed_differences.bed'];
    fi = fopen(fname, 'r');
    Asian_African_mutations = textscan(fi, '%s\t%d\t%d');
    fclose (fi);
    p_CpGs = RoAM_SAM2CpGinfo (filename, Library,chr_lengths, GenomeAssembly, Organism_kind, Organism_spec, 'CHRs', 1:27,'Asian_African_mutations',Asian_African_mutations,'tot_CHRs',CHRs); 
 
else
    p_CpGs = RoAM_SAM2CpGinfo (filename, Library,chr_lengths, GenomeAssembly, Organism_kind, Organism_spec,'TrimEnds',1,'chrFullOrNum','Num'); 
end

as = amSample;
as = as.bundle(p_CpGs);
as = as.set('name',sprintf('I%s',sample),'species','Homo Sapiens',...
        'chr_names',chr_names,'abbrev',sample,'library',Library,...
        'reference',reference);

save(sprintf('u_%s',sample),'as','-v7.3')

%% Diagnose
load (['u_',sample])

HomCov = load('data_coverage', ['HomCov_',sample]);
HomCov = HomCov.(['HomCov_',sample]);
max_coverage = 100;
max_tOct = 0.25;
max_aOga=2*(1/(HomCov/2));

as = as.diagnose('max_tOct',max_tOct,'max_aOga',max_aOga,...
        'max_coverage',max_coverage,'strict',false);

save(sprintf('u_%s',sample),'as','-v7.3')
 
%% Filter and merge
load (['u_',sample])
fas = as.filter();
save(sprintf('f_%s',SAMPLE),'fas');

%% Calculate deamination rate
load (['f_',sample])
if strcmp(Organism_Kind,'Homo Sapiens')
    load 'RefBone5'
end

fas = fas.estimatedrate('method','reference','reference',mms);
save(sprintf('f_%s',SAMPLE),'fas');

%% Reconstruct methylation maps
load (['f_',sample])
fas = fas.reconstructmethylation;
save(sprintf('f_%s',SAMPLE),'fas');


%% Simulate

load (['u_',sample])
load (['f_',sample])

if strcmp(Organism_kind,'Homo Sapiens')
    load 'RefBone5'
end

drate = fas.get('drate');
drate = drate.rate.global;
simu_as = as.simulate(drate,mms,false);
save (sprintf('simu_%s',sample),'simu_as','-v7.3')

HomCov = load('data_coverage', ['HomCov_',sample]);
HomCov = HomCov.(['HomCov_',sample]);
max_aOga=2*(1/(HomCov/2));

simu_as = simu_as.diagnose('max_tOct',max_tOct,'max_aOga',max_aOga,...
'max_coverage',max_coverage,'strict',false);
simu_fas = simu_as.filter();
simu_fas = simu_fas.estimatedrate('method','reference','reference',mms);
simu_fas = simu_fas.reconstructmethylation;
save(sprintf('simu_%s',sample),'simu_fas','-append');

%% All together

HomCov = load('data_coverage', ['HomCov_',sample]);
HomCov = HomCov.(['HomCov_',sample]);
max_coverage = 100;
max_tOct = 0.25;
max_aOga=2*(1/(HomCov/2));

if strcmp(Organism_Kind,'Homo Sapiens')
    load 'RefBone5'
end
    
as = as.diagnose('max_tOct',max_tOct,'max_aOga',max_aOga,...
    'max_coverage',max_coverage,'strict',false);
fas = as.filter();
fas = fas.estimatedrate('method','reference','reference',mms);
fas = fas.reconstructmethylation;

drate = fas.get('drate');
drate = drate.rate.global;
simu_as = as.simulate(drate,mms,false);

simu_as = simu_as.diagnose('max_tOct',max_tOct,'max_aOga',max_aOga,...
'max_coverage',max_coverage,'strict',false);
simu_fas = simu_as.filter();
simu_fas = simu_fas.estimatedrate('method','reference','reference',mms);
simu_fas = simu_fas.reconstructmethylation;


save(sprintf('u_%s',sample),'as','-v7.3');
save(sprintf('f_%s',sample),'fas');
save(sprintf('simu_%s',sample),'simu_as','-v7.3')
save(sprintf('simu_%s',sample),'simu_fas','-append');