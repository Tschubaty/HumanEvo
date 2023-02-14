function show(dm)

% SHOW displays class content.
% --------
% show(dm)
% --------
% Description: display DMRs content.
% Input:       {dm} DMRs array.

% © Liran Carmel
% Classification: Display
% Last revision date: 31-Jul-2017

% constant indentations
PRIMARY_IDENT = 4;
SECONDARY_IDENT = 7;
WIDTH = 80;

% display title
disp(' ');
disp([inputname(1),' = '])
disp(' ');

% use different displays when {dm} is empty
if isscalar(dm) && (isempty(dm) || ~dm.no_samples)
    disp('   Empty DMRs object')
else   % non-empty {dm}
    % display general description of {dm}
    base_samples = dm(1).samples;
    no_samples = dm(1).no_samples;
    fprintf(1,'Comparison between %d samples across %d chromosomes\n',...
        no_samples,length(dm));
    % samples
    str = 'SAMPLES: ';
    for ss = 1:(no_samples-1)
        str = sprintf('%s %s-',str,dm(1).samples{ss});
    end
    str = sprintf('%s%s\n',str,dm(1).samples(no_samples));
    disptext(str,WIDTH,PRIMARY_IDENT,SECONDARY_IDENT);
    % total number of DMRs
    str = 'DMRs:';
    str = sprintf('%s ',);
        if novariables(vsm) == 1
            hypertext = 'variable';
        else
            hypertext = 'variables';
        end
        str = sprintf('%sCOMPOSITION: %d <a href="%s">%s</a>',...
            blanks(PRIMARY_IDENT),novariables(vsm),cmd,hypertext);
        cmd = sprintf('matlab: %s.col_sampleset',inputname(1));
        str = sprintf('%s x %d <a href="%s">samples</a>',...
            str,nosamples(vsm),cmd);
        if nogroupings(vsm) == 1
            hypertext = 'grouping';
        else
            hypertext = 'groupings';
        end
        cmd = sprintf('matlab: %s.col_groupings',inputname(1));
        str = sprintf('%s (%d <a href="%s">%s</a>)',...
            str,nogroupings(vsm),cmd,hypertext);
        disp(str);
        % add description
        if ~isempty(get(vsm,'description'))
            disptext(sprintf('DESCRIPTION: %s',get(vsm,'description')),...
                WIDTH,PRIMARY_IDENT,SECONDARY_IDENT);
        end
        % add source
        if ~isempty(get(vsm,'source'))
            disptext(sprintf('SOURCE: %s',get(vsm,'source')),WIDTH,...
                PRIMARY_IDENT,SECONDARY_IDENT);
        end
    end
else    % if a vector of matrices should be displayed
    for ii = 1:length(vsm)
        str = sprintf('(%d)',ii);
        if novariables(vsm(ii))*nosamples(vsm(ii)) == 0   % default instance
            str = sprintf('%s Empty matrix',str);
        else                                 % non-default instance
            % display name of vsmatrix
            str = sprintf('%s %s (vsmatrix)',str,get(vsm(ii),'name'));
            % display rows and columns information
            str = sprintf('%s: %d variables x %d samples',...
                str,novariables(vsm(ii)),nosamples(vsm(ii)));
            if nogroupings(vsm) == 1
                hypertext = 'grouping';
            else
                hypertext = 'groupings';
            end
            str = sprintf('%s (%d %s)',str,nogroupings(vsm(ii)),hypertext);
        end
        disp(str);
    end
end
disp(' ');