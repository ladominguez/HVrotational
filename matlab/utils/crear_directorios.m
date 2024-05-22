function crear_directorios(rutagrab, estac)

    if ~exist(rutagrab,'dir') 
        mkdir(rutagrab) 
    end

    if ~exist([rutagrab,estac],'dir') 
        mkdir([rutagrab,estac])
    end
       
end