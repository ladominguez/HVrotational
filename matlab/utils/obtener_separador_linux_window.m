function sep = obtener_separador_linux_window()
    if ismac ||  isunix
        sep = '/';
    elseif ispc
        sep = '\';
    else
        disp('Platform not supported')
    end 
