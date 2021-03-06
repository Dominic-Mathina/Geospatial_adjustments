Photo1=[-13.810 38.525         %%Photo coordinates of photo1
        -12.890 9.233
        -12.770 -16.291
         63.754 49.820
         59.290 13.632
         51.541 -10.573
         32.351 31.818
         25.987 7.012];
 x1=Photo1(:,1);
 y1=Photo1(:,2);

Photo2=[-77.425 39.310
        -73.860 10.007
        -73.594 -15.499
         3.328  50.492
         2.061  14.065
        -3.823  -9.889
        -28.976 32.521
        -28.819 7.704];
 x2=Photo2(:,1);
 y2=Photo2(:,2);
 x=[x1;x2];
 y=[y1;y2];
%%
  %Control Points
  CTRL=[2594.79 695.81 2022.36
        2607.84 206.06 2080.31
        3370.35 816.31 1986.42
        3263.70 245.96 2129.39];
    Xc=CTRL(:,1);
    Yc=CTRL(:,2);
    Zc=CTRL(:,3);
    
     %FOCAL LENGTHS1
   f1=165.89;
   f2=165.77;
   f=[f1;f2];
 %%
   %Initial approximate values for photo1
   Xo1=2480;
   Yo1=378;
   Zo1=559;
   om1=0;
   phi1=0;
   K1=0;
   %Initial approximate values for photo2

   Xo2=3033;
   Yo2=377;
   Zo2=469;
   
   om2=0;
   phi2=0;
   K2=0;
   %%%%concatanate these values%%%%
   om=[om1;om2];
   phi=[phi1;phi2];
   K=[K2;K1];
   %%%%%%%%%%%%%%%%
   Xo=[Xo1;Xo2];
   Yo=[Yo1;Yo2];
   Zo=[Zo1;Zo2];
   
   %Initial approximate values for new points 2,5,7 and 8
   Points=[2606 439 2088
           3320 488 2107
           3015 645 2008
           3013 424 2188];
       Xp=Points(:,1);
       Yp=Points(:,1);
       Zp=Points(:,1);
       
       X=[Xc;Xp];
       Y=[Yc;Yp];
       Z=[Zc;Zp];

T=[3;3;9;9;13;13;15;15;19;19;25;25;29;29;31;31];
V=[13;13;14;14;15;15;16;16;17;17;19;19;20;20;21;21;22;22;23;23;25;25];
U=[2;2;5;5;7;7;8;8;10;10;13;13;15;15;16;16];
K=[1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2];
P=[1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2];
for i=1:32
    t=T(i);
i=t;
v=V(i);
j2=v;
u=U(i);
p=P(i);
k=K(i);
    
    
r11=cos(phi(k))*cos(K(k));
 r12=((cos(om(k))*sin(K(k)))+(sin(om(k))*sin(phi(k))*cos(K(k))));
 r13=(sin(om(k))*sin(K(k)))-(cos(om(k))*sin(phi(k))*cos(K(k)));
 r21=-cos(phi(k))*sin(K(k));
 r22=(cos(om(k))*cos(K(k)))-(sin(om(k))*sin(phi(k))*sin(K(k)));
 r23=(sin(om(k))*cos(K(k)))+(cos(om(k))*sin(phi(k))*sin(K(k)));
 r31=sin(phi(k));
 r32=-sin(om(k))*cos(phi(k));
 r33=cos(om(k))*cos(phi(k));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A2(i,0+j2)=x(u)*r31+f(p)*r11;
A2(i+1,0+j2)=y(u)*r31+f(p)*r21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Ya differentials%%
A2(i,1+j2)=-x(u)*r33-f(p)*r13;
A2(i+1,1+j2)-y(u)*r33-f(p)*r23;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Za differentials %%
A2(i,2+j2)=x(u)*r32+f(p)*r12;
A2(i+1,2+j2)=y(u)*r32+f(p)*r22;

end