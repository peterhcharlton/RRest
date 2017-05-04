%% save_or_append_data - save data to either a new file or appended old file
% determines which way to save the data, and saves it.

if exist(savepath, 'file')
    save(savepath, save_name, '-append')
else
    save(savepath, save_name)
end