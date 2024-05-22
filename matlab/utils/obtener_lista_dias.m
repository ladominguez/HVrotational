function listdias = obtener_lista_dias(listreg)

    listdias00 = {};
    for i = 1:length(listreg)
        listdias00{i,1} = listreg{i}(1:8);
     end
     listdias = unique(listdias00);
