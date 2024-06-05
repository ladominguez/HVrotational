function [EWv, NSv, VEv] = remover_media_taper(EWv, NSv, VEv, ptosvent, factap)
    
    EWv = EWv-mean(EWv);
    NSv = NSv-mean(NSv);
    VEv = VEv-mean(VEv);
    tap = repmat(tukeywin(ptosvent,factap),1,length(NSv(1,:)));
    EWv = EWv.*tap;
    NSv = NSv.*tap;
    VEv = VEv.*tap;