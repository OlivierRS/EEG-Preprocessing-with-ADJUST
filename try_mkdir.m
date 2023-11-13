function try_mkdir(folder)
    [folder,~,~] = fileparts(folder);
    if ~isfolder(folder)
        mkdir(folder);
    end
end


