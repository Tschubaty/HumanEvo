function fh = plot(obj,chr,idx_DMR,gc,samples,varargin)

% PLOT plots methylation and supplementary information of a DMR.
% --------------------------------------------------------------
% fh = obj.plot(chr, idx_DMR, gc, samples, i_name, i_value, ...)
% --------------------------------------------------------------
% Description:  plots methylation levels of samples in a certain DMR, as
%               well as additional supplementary information such as gene
%               bodies, CGIs and more.
% Input:        {chr} chromosome where the DMR is located.
%               {idx_DMR} DMR number.
%               {gc} genomic coordinates of CpG positions (GCOORDINATES
%                   object).
%               {samples} cell array of samples. Each element is either an
%                   ancient sample (of type AMSAMPLE) or a modern sample
%                   (of type MMSAMPLE).
%               <{i_name, i_value}> pairs that modify the algorithm.
%                   Currently available:
%                   'orderby' how to order the samples. Can be 'none'
%                       (original order as in {obj}) or 'groups' (def:
%                       'groups').
%                   'genes' GINTERVALS object with gene bodies.
%                   'CGIs' GINTERVAL object with CGIs.
%                   'widenby' number of bases by which to widen the regions
%                       to plot. This number is added to both sides of the
%                       region (def: 0).
% Output:       {fh} figure handle.

% © Liran Carmel
% Classification: Properties
% Last revision date: 14-Dec-2018

% define region
beg = obj.cDMRs(chr).gen_begin(idx_DMR);
fin = obj.cDMRs(chr).gen_end(idx_DMR);
reg = struct('chromosome',obj.chromosomes{chr},'beg',beg,'end',fin);

% identify order
orderby = 'groups';
for ii = 6:2:nargin
    switch str2keyword(varargin{ii-5},6) %#ok<*AGROW>
        case 'orderb'   % instruction: orderby
            orderby = varargin{ii-4};
    end
end

% reorder samples
if strcmp(orderby,'groups')
    idx = [];
    gid = obj.groups.gcn2gid(1);
    gid = gid{1};
    for gg = 1:obj.groups.nogroups
        idx = [idx obj.groups.grp2samp(gid(gg))];
    end
    samples = samples(idx);
end

% call PLOTREGION engine
fh = plotregion(reg,gc,samples,varargin{:});

% add line separating the groups
if strcmp(orderby,'groups')
    hold on
    gs = cumsum(obj.groups.groupsize) + 1;
    for gg = 1:length(gs)-1
        line(xlim,[gs(gg) gs(gg)],'color','k','linewidth',2);
    end
    hold off
end