%% Check_exists - checks to see if the specified field exists in the specified path
% save_name - specified file
% savepath - specified path

if exist(savepath, 'file')
    filecontents = whos('-file', savepath);
    var_names = extractfield(filecontents, 'name');
    if sum(strcmp(var_names, save_name))
        continue
    end
end