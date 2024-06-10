function col = get_colors(itertot)

figure(100)
p = load('progresivo');
p = p.mycmap;
p = colormap(jet);
close(100)

Ncomb = itertot;
porc = length(p(:,1))/Ncomb;
for i = 1:Ncomb
    col(i,:) = p(ceil(porc*i),:);
end
col = flipud(col);
