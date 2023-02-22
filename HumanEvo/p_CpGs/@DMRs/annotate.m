function obj = annotate(obj,gi_CGIs)

% anotate retrieves important data about each DMR
% --------------------------------
% dm = dm.annotate()
% --------------------------------------------------------------------
% Description:  retrieves important data about each DMR.
% Input:        {dm} DMRs object.
% Output:       {dm} a modified cDMRs object.

% © Yoav Mathov & David Gokhman & Liran Carmel
% Classification: RoAM
% Last revision date: 16-Jan-2019

% initialize
ann = [];   % annotation structure

% loop on chromosomes
no_DMRs = obj.noDMRs;
for chr = 1:obj.no_chromosomes
    % return if the chromosome contains no DMRs
    if ~no_DMRs(chr)
        return;
    end
    % initialize annotation structure
    emptycell = cell(no_DMRs(chr),1);
    cann = struct('in_CGI',emptycell);
    % get DMR regions
    regions = obj.cDMRs(chr).getregion();
    % loop on DMRs
    for dd = 1:no_DMRs
        gi = gi_CGIs.overlap(regions(dd));
    end
end


%load island and genebody data
islands = cell (1,length(CHRs));
genebody = cell (1,length(CHRs));
for chr = CHRs
    for dd = 1:obj(chr).no_DMRs
        islands{chr}(dd) = nanmax(TA_ratio_island{chr}((obj(chr).CpG_begin(dd)):(obj(chr).CpG_begin(dd))));
        if nanmax(TA_ratio_genebody{chr}((obj(chr).CpG_begin(dd)):obj(chr).CpG_begin(dd))) == 1
            genebody{chr}(dd) = 1;
        elseif nanmin(TA_ratio_genebody{chr}((obj(chr).CpG_begin(dd)):obj(chr).CpG_begin(dd))) == -1
            genebody{chr}(dd) = -1;
        else
            genebody{chr}(dd) = 0;
        end
    end
end

%####################gene anotation##################################
DMR_list_Annotation = [];
% load UCSC genes
fname = '\\132.64.65.160\David.Gokhman\Human\Methylome\UCSC_Genes.txt';
fid = fopen(fname,'r');
UCSC = textscan(fid,'%d\t%d\t%d\t%s\t%s\t%s %*[^\n]');
fclose(fid);

fprintf(1,'creating promoter table \n')
UCSC_Proms = [];
UCSC_Proms{1} = UCSC{1}; UCSC_Proms{4} = UCSC{4}; UCSC_Proms{5} = UCSC{5}; UCSC_Proms{6} = UCSC{6};
for gene = 1:length(UCSC{1})
    if strcmp(UCSC{6}{gene},'+')
        UCSC_Proms{2}(gene,1) = UCSC{2}(gene)-5000;
        UCSC_Proms{3}(gene,1) = UCSC{2}(gene)+1000;
    elseif strcmp(UCSC{6}{gene},'-')
        UCSC_Proms{2}(gene,1) = UCSC{3}(gene)-1000;
        UCSC_Proms{3}(gene,1) = UCSC{3}(gene)+5000;
    end
end

fprintf(1,'annotating based on %d UCSC genes \n', length(UCSC{1}))

DMR_list_Annotation{1} = []; %all overlapping proms sym
DMR_list_Annotation{2} = []; %all overlapping genes sym

DMR_counter=0;
for chr = CHRs
    for ll = 1:obj(chr).no_DMRs
        DMR_counter = DMR_counter+1;
        AllOverlappingGenes_sym = [];
        AllOverlappingPromoters_sym = [];

        idxOverlap = find(chr == UCSC{1} & ((obj(chr).gen_begin(ll) >= UCSC{2} & obj(chr).gen_begin(ll) <= UCSC{3})  |  (obj(chr).gen_end(ll)+1 >= UCSC{2} & obj(chr).gen_end(ll)+1 <= UCSC{3}) | (obj(chr).gen_begin(ll) <= UCSC{2} & obj(chr).gen_end(ll)+1 >= UCSC{3})));
        idxProms = find(chr == UCSC_Proms{1} & ((obj(chr).gen_begin(ll) >= UCSC_Proms{2} & obj(chr).gen_begin(ll) <= UCSC_Proms{3})  |  (obj(chr).gen_end(ll)+1 >= UCSC_Proms{2} & obj(chr).gen_end(ll)+1 <= UCSC_Proms{3}) | (obj(chr).gen_begin(ll) <= UCSC_Proms{2} & obj(chr).gen_end(ll)+1 >= UCSC_Proms{3})));

        if ~isempty(idxOverlap)
            for ii = 1:length(idxOverlap)
                CurrIdx = idxOverlap(ii);
                AllOverlappingGenes_sym = strcat(AllOverlappingGenes_sym,',',UCSC{4}(CurrIdx));
            end
        end
        if ~isempty(idxProms)
            for ii = 1:length(idxProms)
                CurrIdx = idxProms(ii);
                AllOverlappingPromoters_sym = strcat(AllOverlappingPromoters_sym,',',UCSC_Proms{4}(CurrIdx));
            end
        end

        %remove duplicates
        if ~isempty(AllOverlappingPromoters_sym)
            AllOverlappingPromoters_sym = removeDupsString(AllOverlappingPromoters_sym{1},',');
        end
        DMR_list_Annotation{1}{DMR_counter} = AllOverlappingPromoters_sym;

        if ~isempty(AllOverlappingGenes_sym)
            AllOverlappingGenes_sym = removeDupsString(AllOverlappingGenes_sym{1},',');
        end
        DMR_list_Annotation{2}{DMR_counter} = AllOverlappingGenes_sym;
    end
end


DMR_list_Annotation{3} = []; %closest dist upstream
DMR_list_Annotation{4} = []; %closest gene upstream
DMR_list_Annotation{5} = []; %closest dist downstream
DMR_list_Annotation{6} = []; %closest gene downstream
DMR_list_Annotation{7} = []; %mindist to TSS
DMR_list_Annotation{8} = {}; %closest gene to TSS


%closest genes
DMR_counter=0;
fprintf(1,'annotating closest gene \n')

for chr = CHRs
    for ll = 1:obj(chr).no_DMRs
        DMR_counter = DMR_counter+1;

        ClosestDistUpstream = -999999999;
        ClosestDistDownstream = 999999999;
        idxUp = 0;
        idxDown = 0;
        for ii = 1:length(UCSC{1})
            UCSCchr = UCSC{1}(ii); %UCSC_TSS_raw{ii,1}
            if UCSCchr ~= chr
                continue
            else
                if strcmp(UCSC{6}(ii),'+')
                    UCSCTSS = UCSC{2}(ii);
                else
                    UCSCTSS = UCSC{3}(ii);
                end
                dist1 = obj(chr).gen_begin(ll) - UCSCTSS;
                dist2 = obj(chr).gen_end(ll)+1 - UCSCTSS;

                if abs(dist1) < abs(dist2)
                    dist = dist1;
                else
                    dist = dist2;
                end

                UCSCStrand = UCSC{6}(ii);
                if  obj(chr).gen_begin(ll) <= UCSCTSS && obj(chr).gen_end(ll)+1 >= UCSCTSS
                    ClosestDistUpstream = 0;
                    ClosestDistDownstream = 0;
                    idxUp = ii;
                    idxDown = ii;
                elseif (dist < 0) && dist > ClosestDistUpstream && strcmp(UCSCStrand,'+')
                    ClosestDistUpstream = dist;
                    idxUp = ii;
                elseif (dist > 0) && dist < ClosestDistDownstream && strcmp(UCSCStrand,'+')
                    ClosestDistDownstream = dist;
                    idxDown = ii;
                elseif (dist < 0) && abs(dist) < ClosestDistDownstream && strcmp(UCSCStrand,'-')
                    ClosestDistDownstream = abs(dist);
                    idxDown = ii;
                elseif (dist > 0) && -1*dist > ClosestDistUpstream && strcmp(UCSCStrand,'-')
                    ClosestDistUpstream = -1*dist;
                    idxUp = ii;
                end
            end
        end


        if idxDown == 0 && idxUp == 0 %can't be! added because of a bug I couldn't find
            DMR_list_Annotation{3}(DMR_counter) = ClosestDistUpstream;
            DMR_list_Annotation{4}{DMR_counter} = 'noGeneDownstream';
            DMR_list_Annotation{5}(DMR_counter) = ClosestDistDownstream;
            DMR_list_Annotation{6}{DMR_counter} = 'noGeneDownstream';
        elseif idxDown == 0
            DMR_list_Annotation{3}(DMR_counter) = ClosestDistUpstream;
            DMR_list_Annotation{4}{DMR_counter} = UCSC{4}{idxUp};
            DMR_list_Annotation{5}(DMR_counter) = ClosestDistDownstream;
            DMR_list_Annotation{6}{DMR_counter} = 'noGeneDownstream';
        elseif idxUp == 0
            DMR_list_Annotation{3}(DMR_counter) = ClosestDistUpstream;
            DMR_list_Annotation{4}{DMR_counter} = 'noGeneUpstream';
            DMR_list_Annotation{5}(DMR_counter) = ClosestDistDownstream;
            DMR_list_Annotation{6}{DMR_counter} = UCSC{4}{idxDown};

        else
            DMR_list_Annotation{3}(DMR_counter) = ClosestDistUpstream;
            DMR_list_Annotation{4}{DMR_counter} = UCSC{4}{idxUp};
            DMR_list_Annotation{5}(DMR_counter) = ClosestDistDownstream;
            DMR_list_Annotation{6}{DMR_counter} = UCSC{4}{idxDown};
        end
    end
end

%minDist to TSS and Closest TSS 
AbsDist = abs(DMR_list_Annotation{3}) - abs(DMR_list_Annotation{5});
idxUps = AbsDist <= 0;
idxDow = AbsDist > 0;
DMR_list_Annotation{7}(idxUps) = DMR_list_Annotation{3}(idxUps);
DMR_list_Annotation{7}(idxDow) = DMR_list_Annotation{3}(idxDow);

DMR_list_Annotation{8}(idxUps) = DMR_list_Annotation{4}(idxUps);
DMR_list_Annotation{8}(idxDow) = DMR_list_Annotation{6}(idxDow);



for ii = 1:length(DMR_list_Annotation)
    DMR_list_Annotation{ii} = DMR_list_Annotation{ii}';
end


%###################Calc methylation######################
fprintf(1,'Calculating methylation \n');
DMR_meth = []; 

col_counter = 1;


%General group methylation
if ~isempty  (obj(chr).methylation)
    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;

            DMR_meth(DMR_counter,col_counter) = obj(chr).methylation(1,dd);
            DMR_meth(DMR_counter,col_counter+1) = obj(chr).methylation(2,dd);        
        end
    end
    col_counter = col_counter+2;
end

for smp = 1:length(obj(1).samples)
    load (['meth_ancient_',obj(1).kind],['TA_convperc_',obj(1).samples{smp}])
    temp = eval(['TA_convperc_',obj(1).samples{smp}]);

    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;
            begTA = obj(chr).CpG_begin(dd);
            finTA = obj(chr).CpG_end(dd)+1;
            DMR_meth(DMR_counter,col_counter) = nanmean(temp{chr}(begTA:finTA));

        end
    end

    clearvars(['TA_convperc_',obj(1).samples{smp}])
    col_counter = col_counter+1;
end

if strcmp (obj(1).kind,'Human')   
    %Chimp WGBS bone
    load meth_ChimpBoneWGBS ChimpBoneWGBS
    ChimpBoneWGBS = ChimpBoneWGBS{1};
    load meth_ChimpBoneWGBS_smoothed ChimpBoneWGBS_Smoothed
    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;
            begTA = obj(chr).CpG_begin(dd);
            finTA = obj(chr).CpG_end(dd)+1;
            DMR_meth(DMR_counter,col_counter) = nanmean(ChimpBoneWGBS{chr}(begTA:finTA));
            DMR_meth(DMR_counter,col_counter+1) = nanmean(ChimpBoneWGBS_Smoothed{1}{chr}(begTA:finTA));
        end
    end
    clearvars ChimpBoneWGBS_Smoothed


    %Chimp RRBS bone
    load meth_ChimpBoneRRBS ChimpBoneRRBS
    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;
            begTA = obj(chr).CpG_begin(dd);
            finTA = obj(chr).CpG_end(dd)+1;

            idxRRBS_Chimp = ~isnan(ChimpBoneRRBS{1}{chr}(begTA:finTA)) | ~isnan(ChimpBoneRRBS{2}{chr}(begTA:finTA));
            posRRBS = begTA:finTA; posRRBS = posRRBS(idxRRBS_Chimp);
            meanRRBS_Chimp = nanmean([ChimpBoneRRBS{1}{chr}(begTA:finTA); ChimpBoneRRBS{2}{chr}(begTA:finTA)]);
            ratioWGBS2RRBS = nanmean(ChimpBoneWGBS{chr}(begTA:finTA))/nanmean(ChimpBoneWGBS{chr}(posRRBS)); %nanmean(ChimpBoneWGBS{chr}(posRRBS)) / meanRRBS_Chimp;
            if ~isnan(ratioWGBS2RRBS)
                DMR_meth(DMR_counter,col_counter+2) = min([100,ratioWGBS2RRBS * meanRRBS_Chimp]); %chimp RRBS normalized
            else
                DMR_meth(DMR_counter,col_counter+2) = nan;
            end
        end
    end
    clearvars ChimpBoneRRBS

    %Chimp 850K
    load meth_ChimpBone850K

    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;
            begTA = obj(chr).CpG_begin(dd);
            finTA = obj(chr).CpG_end(dd)+1;

            %annoate 850K meth chimp bone
            probe_idx = find(coors_850K_chimp_TA(:,1) == chr & coors_850K_chimp_TA(:,2) >= begTA & coors_850K_chimp_TA(:,2) <= finTA);
            DMR_meth(DMR_counter,col_counter+3) = nan;
            DMR_meth(DMR_counter,col_counter+4) = nan;
            DMR_meth(DMR_counter,col_counter+5) = nan;
            DMR_meth(DMR_counter,col_counter+6) = nan;
            DMR_meth(DMR_counter,col_counter+7) = 0;
            DMR_meth(DMR_counter,col_counter+8) = nan;
            DMR_meth(DMR_counter,col_counter+9) = nan;
            DMR_meth(DMR_counter,col_counter+10) = nan;

            if ~isempty(probe_idx)
                meanMeth850K = nanmean(nanmean(Meth850K_chimp(probe_idx,:) * 100));
                minMeth850K = nanmin(nanmean(Meth850K_chimp(probe_idx,:) * 100,1));
                maxMeth850K = nanmax(nanmean(Meth850K_chimp(probe_idx,:) * 100,1));
                DMR_meth(DMR_counter,col_counter+3) = meanMeth850K;
                DMR_meth(DMR_counter,col_counter+4) = minMeth850K;
                DMR_meth(DMR_counter,col_counter+5) = maxMeth850K;
                DMR_meth(DMR_counter,col_counter+6) = nanstd(nanmean(Meth850K_chimp(probe_idx,:) * 100,1));
                DMR_meth(DMR_counter,col_counter+7) = length(probe_idx);
                coors_probes = coors_850K_chimp(probe_idx,2);
                coors_probes_TA = coors_850K_chimp_TA(probe_idx,2);

                if ~isempty(coors_probes_TA)
                    if ~isnan(nanmean(ChimpBoneWGBS{chr}(coors_probes_TA)))
                        ratioWGBS2850K = nanmean(ChimpBoneWGBS{chr}(begTA:finTA))/nanmean(ChimpBoneWGBS{chr}(coors_probes_TA)); %nanmean(ChimpBoneWGBS{chr}(coors_probes_TA)) / meanMeth850K;
                        DMR_meth(DMR_counter,col_counter+8) = min([100,ratioWGBS2850K * meanMeth850K]); %normalized mean meth
                        DMR_meth(DMR_counter,col_counter+9) = min([100,ratioWGBS2850K * minMeth850K]);
                        DMR_meth(DMR_counter,col_counter+10) = min([100,ratioWGBS2850K * maxMeth850K]);
                    end
                end
            end
        end
    end
    clearvars ChimpBoneWGBS coors_850K_chimp Meth850K_chimp probes_850K_chimp coors_850K_chimp_TA

    %Human WGBS Bone
    MethMat4_smoothed = load('meth_BoneWGBS_smoothed','MethMat4'); MethMat4_smoothed = MethMat4_smoothed.MethMat4;
    MethMat5_smoothed = load('meth_BoneWGBS_smoothed','MethMat5'); MethMat5_smoothed = MethMat5_smoothed.MethMat5;
    load meth_BoneWGBS MethMat4 MethMat5
    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;
            begTA = obj(chr).CpG_begin(dd);
            finTA = obj(chr).CpG_end(dd)+1;

            DMR_meth(DMR_counter,col_counter+11) = nanmean(MethMat4{chr}(begTA:finTA)); %bone meth samp4
            DMR_meth(DMR_counter,col_counter+12) = nanmean(MethMat5{chr}(begTA:finTA)); %bone meth samp5 (the better sample)
            DMR_meth(DMR_counter,col_counter+13) = nanmean(MethMat4_smoothed{chr}(begTA:finTA));
            DMR_meth(DMR_counter,col_counter+14) = nanmean(MethMat5_smoothed{chr}(begTA:finTA));
        end
    end
    clearvars MethMat4_smoothed MethMat5_smoothed



    %RRBS Osteoblast
    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;
            begTA = obj(chr).CpG_begin(dd);
            finTA = obj(chr).CpG_end(dd)+1;

            idxRRBS_Osteo = find(~isnan(TA_ratio_meth_Osteo{chr}(begTA:finTA)));
            posRRBS = [begTA:finTA]; posRRBS = posRRBS(idxRRBS_Osteo);
            meanRRBS_Osteo = nanmean(TA_ratio_meth_Osteo{chr}(begTA:finTA));
            if ~isnan(nanmean([MethMat4{chr}(posRRBS); MethMat5{chr}(posRRBS)]))
                ratioWGBS2RRBS = nanmean([MethMat4{chr}(begTA:finTA); MethMat5{chr}(begTA:finTA)]) / nanmean([MethMat4{chr}(posRRBS); MethMat5{chr}(posRRBS)]);
                DMR_meth(DMR_counter,col_counter+15) = min([100,ratioWGBS2RRBS * meanRRBS_Osteo]);
            else
                DMR_meth(DMR_counter,col_counter+15) = nan;
            end
        end
    end

    %450K Human
    load meth_bone_450K Meth450K coors_450K
    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;
            begTA = obj(chr).CpG_begin(dd);
            finTA = obj(chr).CpG_end(dd)+1;

            probe_idx = find(coors_450K(:,1) == chr & coors_450K(:,2) >= begTA & coors_450K(:,2) <= finTA);
            DMR_meth(DMR_counter,col_counter+16) = nan;
            DMR_meth(DMR_counter,col_counter+17) = nan;
            DMR_meth(DMR_counter,col_counter+18) = nan;
            DMR_meth(DMR_counter,col_counter+19) = nan;
            DMR_meth(DMR_counter,col_counter+20) = 0;
            DMR_meth(DMR_counter,col_counter+21) = nan;
            DMR_meth(DMR_counter,col_counter+22) = nan;
            DMR_meth(DMR_counter,col_counter+23) = nan;

            if ~isempty(probe_idx)
                meanMeth450K = nanmean(nanmean(Meth450K(probe_idx,:) * 100));
                minMeth450K = nanmin(nanmean(Meth450K(probe_idx,:) * 100,1));
                maxMeth450K = nanmax(nanmean(Meth450K(probe_idx,:) * 100,1));
                DMR_meth(DMR_counter,col_counter+16) = meanMeth450K;
                DMR_meth(DMR_counter,col_counter+17) = minMeth450K;
                DMR_meth(DMR_counter,col_counter+18) = maxMeth450K;
                DMR_meth(DMR_counter,col_counter+19) = nanstd(nanmean(Meth450K(probe_idx,:) * 100,1));
                DMR_meth(DMR_counter,col_counter+20) = length(probe_idx);
                coors_probes = coors_450K(probe_idx,2);
                coors_probes_TA = [];
                for cc = 1:length(coors_probes)
                    idx = find(TA_ratio_coor{chr} == coors_probes(cc));
                    coors_probes_TA = [coors_probes_TA idx];
                end
                if ~isempty(coors_probes_TA) && ~isnan(nanmean([nanmean(MethMat5{chr}(coors_probes_TA)) nanmean(MethMat4{chr}(coors_probes_TA))]))
                    ratioWGBS2450K = nanmean([nanmean(MethMat5{chr}(begTA:finTA)) nanmean(MethMat4{chr}(begTA:finTA))])/nanmean([nanmean(MethMat5{chr}(coors_probes_TA)) nanmean(MethMat4{chr}(coors_probes_TA))]);
                    DMR_meth(DMR_counter,col_counter+21) = min([100,ratioWGBS2450K * meanMeth450K]); %normalized mean meth
                    DMR_meth(DMR_counter,col_counter+22) = min([100,ratioWGBS2450K * minMeth450K]);
                    DMR_meth(DMR_counter,col_counter+23) = min([100,ratioWGBS2450K * maxMeth450K]);
                end
            end
        end
    end
    clearvars MethMat5 MethMat4  Meth450K coors_450K


    %Human and Chimp brains (+24 and +25 will be the averages        
    load meth_HumanBrain Unmerged_HumanBrain
    load meth_ChimpBrain Unmerged_ChimpBrain
    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;
            begTA = obj(chr).CpG_begin(dd);
            finTA = obj(chr).CpG_end(dd)+1;

            DMR_meth(DMR_counter,col_counter+26) = nanmean(Unmerged_HumanBrain{1}{chr}(begTA:finTA));
            DMR_meth(DMR_counter,col_counter+27) = nanmean(Unmerged_HumanBrain{2}{chr}(begTA:finTA));
            DMR_meth(DMR_counter,col_counter+28) = nanmean(Unmerged_HumanBrain{3}{chr}(begTA:finTA));
            DMR_meth(DMR_counter,col_counter+29) = nanmean(Unmerged_ChimpBrain{1}{chr}(begTA:finTA));
            DMR_meth(DMR_counter,col_counter+30) = nanmean(Unmerged_ChimpBrain{2}{chr}(begTA:finTA));
            DMR_meth(DMR_counter,col_counter+31) = nanmean(Unmerged_ChimpBrain{3}{chr}(begTA:finTA));

            DMR_meth(DMR_counter,col_counter+24) = nanmean(DMR_meth(DMR_counter,col_counter+26:col_counter+28));
            DMR_meth(DMR_counter,col_counter+25) = nanmean(DMR_meth(DMR_counter,col_counter+29:col_counter+31));
        end
    end
    clearvars Unmerged_HumanBrain Unmerged_ChimpBrain


    %Blood
    load meth_ApeBloodWGBS MethMat_chimp MethMat_gorilla MethMat_human MethMat_orangutan
    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;
            begTA = obj(chr).CpG_begin(dd);
            finTA = obj(chr).CpG_end(dd)+1;        

            DMR_meth(DMR_counter,col_counter+32) = nanmean(MethMat_human{chr}(begTA:finTA));
            DMR_meth(DMR_counter,col_counter+33) = nanmean(MethMat_chimp{chr}(begTA:finTA));
            DMR_meth(DMR_counter,col_counter+34) = nanmean(MethMat_gorilla{chr}(begTA:finTA));
            DMR_meth(DMR_counter,col_counter+35) = nanmean(MethMat_orangutan{chr}(begTA:finTA));
        end
    end
    clearvars MethMat_chimp MethMat_gorilla MethMat_human MethMat_orangutan


    %Meissner and Roadmap tissues
    load meth_Meissner
    load meth_Roadmap TA_ratio_Roadmap_primary TA_ratio_Roadmap_primary_names
    DMR_counter = 0;
    for chr = CHRs
        for dd = 1:obj(chr).no_DMRs 
            DMR_counter=DMR_counter+1;
            begTA = obj(chr).CpG_begin(dd);
            finTA = obj(chr).CpG_end(dd)+1;   


            methTemp = [];
            for tissue = 1:length(TA_ratio_Roadmap_primary)
                methTemp(tissue) = nanmean(TA_ratio_Roadmap_primary{tissue}{chr}(begTA:finTA));
            end
            for tissue = length(TA_ratio_Roadmap_primary)+1:length(TA_ratio_Roadmap_primary)+length(MeissnerFull)
                methTemp(tissue) = nanmean(MeissnerFull{tissue - length(TA_ratio_Roadmap_primary)}{chr}(begTA:finTA));
            end
            RoadMies{1}(DMR_counter) = nanmean(methTemp);
            RoadMies{2}(DMR_counter) = nanstd(methTemp);
            RoadMies{3}(DMR_counter) = nanmin(methTemp);
            idxMin = find(methTemp == nanmin(methTemp));
            if idxMin(1)<=length(TA_ratio_Roadmap_primary)
               RoadMies{4}(DMR_counter) = TA_ratio_Roadmap_primary_names(idxMin(1));
            else
               RoadMies{4}(DMR_counter) = Names_MeissnerFull(idxMin(1)-length(TA_ratio_Roadmap_primary));
            end
            RoadMies{5}(DMR_counter) = nanmax(methTemp);
            idxMax = find(methTemp == nanmax(methTemp));
            if idxMax(1)<=length(TA_ratio_Roadmap_primary)
                RoadMies{6}(DMR_counter) = TA_ratio_Roadmap_primary_names(idxMax(1));
            else
                RoadMies{6}(DMR_counter) = Names_MeissnerFull(idxMax(1)-length(TA_ratio_Roadmap_primary));
            end
        end
    end     
end


%##########################     Write file      ##############
fprintf(1,'Writing file \n');

if strcmp (obj(1).kind,'Human')
    titles_beg = 'chr\tgenomic start\tgenomic end\timax\tlength\tCpG start\tCpG end\tNo. CpG\tCGI\tgenebody\tall promoter genes\tall ovelapping genes\tupDist to TSS\tgene\tdownDist to TSS\tgene\tminDist\tclosest gene\tgroup1 meth\tgroup2 meth\t';
    regexp_beg = '%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%f\t%f\t';
    titles_fin = 'Chimp meth\tChimp meth smoothed\tChimp RRBS\tmean 850K chimp\tmin 850K chimp\tmax 850K chimp\tSTD 850K chimp\t#pos 850K chimp\tmean 850K chimp norm\tmin 850K chimp norm\tmax 850K chimp norm\tBone meth4\tBone meth5\tBone meth4_smoothed\tBone meth5_smoothed\tOsteoblast RRBS\tmean 450K\tmin 450K\tmax 450K\tSTD 450K\t#pos 450K\tmean 450K norm\tmin 450K norm\tmax 450K norm\tMeanHumanBrain\tMeanChimpBrain\tHumanBrain1\tHumanBrain2\tHumanBrain3\tChimpBrain1\tChimpBrain2\tChimpBrain3\tHumanBlood\tChimpBlood\tGorillaBlood\tOrangutanBlood\tMeanofTissues\tSTDofTissues\tmin value\tmin tissue\tmax value\tmax tissues\n';
    regexp_fin = '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%s\n';
    titles_meth = [];
    regexp_meth = [];
    for tt = 1:length (obj(1).samples)
        titles_meth = strcat(titles_meth,obj(1).samples{tt}, ' meth\t');
        regexp_meth = strcat(regexp_meth,'%d\t');
    end
    titles = strcat(titles_beg, titles_meth, titles_fin);
    reg = [regexp_beg, regexp_meth, regexp_fin];

    filename = '\\132.64.65.160\David.Gokhman\AnnoTable_temp.bed'; 
    fid = fopen(filename,'w');
    fprintf(fid,titles);

    DMR_counter = 0;
    for chr = CHRs
        for ll = 1:obj(chr).no_DMRs
            DMR_counter = DMR_counter+1;
            fprintf(fid,reg,chr,obj(chr).gen_begin(ll),obj(chr).gen_end(ll),...
            obj(chr).max_Qt(ll),obj(chr).no_bases(ll),obj(chr).CpG_begin(ll),...
            obj(chr).CpG_end(ll), obj(chr).no_CpGs(ll),islands{chr}(ll),genebody{chr}(ll),...
            DMR_list_Annotation{1}{DMR_counter},DMR_list_Annotation{2}{DMR_counter},...
            DMR_list_Annotation{3}(DMR_counter),DMR_list_Annotation{4}{DMR_counter},...
            DMR_list_Annotation{5}(DMR_counter),DMR_list_Annotation{6}{DMR_counter},...
            DMR_list_Annotation{7}(DMR_counter),DMR_list_Annotation{8}{DMR_counter},...
            DMR_meth(DMR_counter,1:end),...
            RoadMies{1}(DMR_counter),RoadMies{2}(DMR_counter),RoadMies{3}(DMR_counter),...
            RoadMies{4}{DMR_counter},RoadMies{5}(DMR_counter),RoadMies{6}{DMR_counter});
        end
    end
    fclose(fid);
    fprintf(1, 'The output was saved as AnnoTable_temp.txt \n');
end


end
