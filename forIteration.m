for i=1:99;
    for i=1:33;
    i=ee(i);
    f=ff(i);
    Y(i,1)=dy1(f);
    end
    for i=1:33;
    i=ee(i);
    f=ff(i);
    X(i,1)=dx1(f);
    end
    for i=1:33'
    i=ee(i);
    f=ff(i);
    Z(i,1)=dz1(f);
    end
end
Yobsv=Y;
Xobsv=X;
Zobsv=Z;