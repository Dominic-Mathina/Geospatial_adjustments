A = [1100.64 1431.09 1 0;
    1431.09 -1100.64 0 1;
    1678.39 254.156 1 0;
    254.15 -1678.39 0 1];

L = [632.17;121.45;355.2;-642.07];

d = inv(A)

soln = d*L

Parameters = [soln(1,:) soln(2,:);
              -soln(2,1) soln(1,:)]

Translations = [soln(3,:);soln(4,:)]

B = [1304.81; 596.37];

%%%%%%%%%Solution Men %%%%%%%%

s1 = B-Translations;
Param = inv(Parameters);

Coord = Param*s1