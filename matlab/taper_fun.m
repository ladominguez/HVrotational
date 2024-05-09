function [Xtap] = taper_fun(x,perce,ini,fin)
Ntot = length(x);
npt = round(Ntot*(perce/100.));
if ini == 1
    for i = 1:npt
        dth = (npt-i)*pi/(npt-1);
        a = 0.5*(1.0+cos(dth));
        x(i) = x(i)*a;
    end
end
if fin == 1
    for i = 0:(npt-1)
        dth = i*pi/(npt-1);
        a = 0.5*(1.0+cos(dth));
        bb = Ntot-(npt-1)+i;
        x(bb) = x(bb)*a; 
    end
end
Xtap = x;
