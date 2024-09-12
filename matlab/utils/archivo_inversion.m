function [HVesc,fsesc] = archivo_inversion(HV,f,f1,f2,paso,fs)

HVmeansmooth = suavfrec(HV,f,fs);
i1 = find(f>=f1,1,'first');
i2 = find(f>=f2,1,'first');
HVesc = HVmeansmooth(i1:paso:i2);
fsesc = f(i1:paso:i2);
figure(400)
semilogx(f,HV,'b','linewidth',2); hold on
semilogx(f,HVmeansmooth,'g','linewidth',2); hold on
semilogx(fsesc,HVesc,'r','linewidth',2); hold on
% saveas(gcf,[rutagrab,separador,'HV_',estac,'.png']);
% close(400)
