function rel_els = find_rel_subjs(up, group_type)
%FIND_REL_SUBJS identifies the subjects for each group (YHVs, EHVs and
%OHVs)
%
%	Inputs:
%		db_data         demographic data extracted from the Vortal database
%       up              universal parameters structure
%
%	Outputs:
%       rel_els         a structure of subject els for each group
%

if strcmp(group_type, 'group')
    groups = up.paramSet.groups;
    group_names = unique(groups);
elseif strcmp(group_type, 'sub_group')
    groups = up.paramSet.sub_groups;
    group_names = unique(groups);
end

for group_no = 1 : length(group_names)
    eval(['rel_els.' group_names{group_no} ' = find(strcmp(groups, group_names{group_no} ));']);    
end

%% Old Code
% 
% loaded_info = load(up.paths.db_data);
% 
% temp.group = loaded_info.db_data.group(loaded_info.rel_els);
% 
% if length(unique(temp.group)) == 2
%     rel_els.yhv = find(temp.group == 1);
%     rel_els.ehv = find(temp.group == 2);
%     rel_els.global = 1:length(temp.group);
% elseif length(unique(temp.group)) == 3
%     rel_els.bw = find(temp.group == 1);
%     rel_els.am = find(temp.group == 2);
%     rel_els.fm = find(temp.group == 3);
%     rel_els.global = 1:length(temp.group);
% else
%     fprint('Unexpected Number of Groups');
% end

end
