function [Nv, EWv, NSv, VEv, fechahmsvent]=ventanas_efectivas(EWv, NSv, VEv, ventok, fechahmsvent, Nvent)

    EWv = EWv(:,ventok);
    NSv = NSv(:,ventok);
    VEv = VEv(:,ventok);
    fechahmsvent = fechahmsvent(ventok);
    
    ventNOefectiv = unique([find(sum(abs(EWv))==0) find(sum(abs(NSv))==0) find(sum(abs(VEv))==0)]);
    ventefectiv = 1:Nvent;
    ventefectiv(ventNOefectiv) = [];
    Nv = length(ventefectiv);
    EWv = EWv(:,ventefectiv);
    NSv = NSv(:,ventefectiv);
    VEv = VEv(:,ventefectiv);
    fechahmsvent = fechahmsvent(ventefectiv);