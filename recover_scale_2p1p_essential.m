function s=recover_scale_2p1p_essential(Rji,tji_norm,Tj2j,xi2,xj2)
    % compute s using Rji tji_norm and xiq xj1 Tjj
    % so tji=s*tji_norm
    Rj2j=Tj2j(1:3,1:3);
    tj2j=Tj2j(1:3,4);
    R2=Rj2j*Rji;
    t2A=Rj2j*tji_norm;
    t2b=tj2j;
    r1=mean([R2(1,1),R2(3,3)]);
    r2=mean([R2(1,3),-R2(3,1)]);
    x1=xi2(1);y1=xi2(2);
    x2=xj2(1);y2=xj2(2);
    A1=r2*x1*y2+y1-r1*y2;
    A3=r1*x1*y2-x2*y1+r2*y2;
    a1=t2A(1);
    a3=t2A(3);
    b1=t2b(1);
    b3=t2b(3);
    s=-(A1*b1+A3*b3)/(A1*a1+A3*a3);
end