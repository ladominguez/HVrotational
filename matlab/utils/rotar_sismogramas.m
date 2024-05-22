function [EW, NS, VE, teta] = rotar_sismogramas(ESTR, tetarot, Nteta)

teta = tetarot(Nteta);
EW = ESTR.EW*cosd(teta)+ESTR.NS*sind(teta);  %longitudinal
NS =-ESTR.EW*sind(teta)+ESTR.NS*cosd(teta);  %transversal
VE = ESTR.VE;