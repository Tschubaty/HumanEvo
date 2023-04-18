
% v 01 04.04.2023
% read methylation files
%
%
%% load liran meta data
%meta = readtable('Meth samples - human, apes, mammoth.xlsx','Sheet');
meta = readtable('Meth samples - human, apes, mammoth.xlsx','Sheet','Ancient', ...
    'Format', "auto");
%'%s %s %s %f %f %s %s %s %s %s %s %s %s %s %s %s %f %f %d %s %s %s');



% delete empty rows
empty_rows = cellfun('isempty', meta.Sample);
meta=meta(~empty_rows ,:);

%% load Allen Ancient Genome Diversity Project from https://reich.hms.harvard.edu/ancient-genome-diversity-project
AGDP_meta = readtable('AGDP.metadata.xlsx');

%% availble data riech input

available_data = {};

name_list = AGDP_meta.I_ID;

fprintf('check %d samples of Reich if found in data folder\n',numel(name_list))

for i=1:numel(name_list )

    pattern =  strcat('p_CpGs/f_',erase(name_list{i},lettersPattern),'.mat');
    % strcat('p_CpGs/*',erase(name_list {i},lettersPattern),'*');
    d=dir(pattern);


    % print out
    %         formatSpec = '%d matches found for %s \n';
    %         fprintf(formatSpec,numel(d),name_list {i});

    if numel(d) > 0
        available_data{end +1} = name_list {i};

        for j=1:numel(d)
            fname=d(j).name;
            % process each file here...
            % disp(fname);
        end

    else

        %fprintf("search pattern: %s \n",pattern);
    end


end


fprintf(' %d out of %d samples of Reich found in data folder\n',numel(available_data),numel(name_list))

%% availble data liarn input

available_data_liran_tabel = {};

name_list = meta.Code;

fprintf('check %d samples of liran if found in data folder\n',numel(name_list))

for i=1:numel(name_list )

    pattern =  strcat('p_CpGs/*',name_list{i},'*');
    %pattern =  strcat('p_CpGs/*',erase(name_list{i},"I"),'*');
    % strcat('p_CpGs/*',erase(name_list {i},lettersPattern),'*');
    d=dir(pattern);


    % print out

    %formatSpec = '%d matches found for %s \n';
    %fprintf(formatSpec,numel(d),name_list {i});

    if numel(d) > 0 &&  ~isempty(name_list{i})
        available_data_liran_tabel{end +1} = name_list {i};

        for j=1:numel(d)
            fname=d(j).name;
            % process each file here...
            disp(fname);
        end
    else

        fprintf("search pattern: %s \n",pattern);
    end


end

fprintf(' %d out of %d samples of liran  found in data folder\n',numel(available_data_liran_tabel),numel(name_list))


%% check david reichs data sets

found_in_both = {};
not_found_in_liran = {};

for i=1:numel(available_data)
    index  = contains(available_data_liran_tabel,erase(available_data{i},lettersPattern));

    if any(index)
        found_in_both{end +1} = available_data{i};
    else
         not_found_in_liran{end +1} = available_data{i};
    end

end

fprintf(' %d out of %d samples of available reich sets  found in data folder\n',numel(found_in_both ),numel(available_data))

%% check liran sets 

not_found_in_reich = {};

for i=1:numel(available_data_liran_tabel)
    index  = contains(available_data,erase(available_data_liran_tabel{i},lettersPattern));

    if ~any(index)
        not_found_in_reich{end +1} = available_data_liran_tabel{i};
    end

end

fprintf(' %d out of %d samples of available liran sets  found in reichs website \n',numel(found_in_both ),numel(available_data_liran_tabel))



