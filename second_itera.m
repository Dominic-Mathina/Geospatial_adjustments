clear all;
format long;
%%%% %%%%%%%%%%%%% Authors-cum-Engineers %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MWONGELA D MATHINA: F19/1707/2013
%%% KARERU PAUL       : F19/    /2013


%%
Photo1=[-13.810 38.525         %%Photo coordinates of photo1
        -12.890 9.233
        -12.770 -16.291
         63.754 49.820
         59.290 13.632
         51.541 -10.573
         32.351 31.818
         25.987 7.012];
     P1=[Photo1(1,:) Photo1(2,:) Photo1(3,:) Photo1(4,:) Photo1(5,:) Photo1(6,:) Photo1(7,:) Photo1(8,:)];
 x1=Photo1(:,1);
 y1=Photo1(:,2);
  L1=P1';%%%%%combine x and y respe %%

Photo2=[-77.425 39.310          %%Photo coordinates of photo2     
        -73.860 10.007
        -73.594 -15.499
         3.328  50.492
         2.061  14.065
        -3.823  -9.889
        -28.976 32.521
        -28.819 7.704];
    
     P2=[Photo2(1,:) Photo2(2,:) Photo2(3,:) Photo2(4,:) Photo2(5,:) Photo2(6,:) Photo2(7,:) Photo2(8,:)];
 x2=Photo2(:,1);
 y2=Photo2(:,2);
   L2=P2'; %%%%%combine x and y respe %%
    LL=[L1;L2];%%%%%%The L matrix %%%%%%%
    L=[-13.81
38.525
-77.425
39.31
-12.89
9.233
-73.86
10.007
-12.77
-16.291
-73.594
-15.499
63.754
49.82
3.328
50.492
59.29
13.632
2.061
14.065
51.541
-10.573
-3.823
-9.889
32.351
31.818
-28.976
32.521
25.987
7.012
-28.819
7.704
];
 x=[x1;x2];
 y=[y1;y2];
%%
 %Control Points
  CTRL=[2594.79 695.81 2022.36
        2607.84 206.06 2080.31
        3310.35 816.31 1986.42
        3263.10 245.96 2129.39];
    Xc=CTRL(:,1);
    Yc=CTRL(:,2);
    Zc=CTRL(:,3);
    
 %%
 
 
 
 %Initial approximate values for photo1
 
   Xo1=2.480676636448541*1000;
   Yo1=0.376727377616105*1000;
   Zo1=0.560306913739147*1000;
   om1= -0.000016003109288*1000;
   phi1=-0.000160569580307*1000;
   K1= -0.000004753020711*1000;
   %Initial approximate values for photo2

   Xo2= 3.035523460007414*1000;
   Yo2=0.382533721682260*1000;
   Zo2= 0.464637039111426*1000;
   om2=-0.000019397086137*1000;
   phi2=-0.000159150155937*1000;
   K2=-0.000001266625458*1000;
   %%%%%%%%%%concatanate these values%%%%%%%%%%%%
   om=[om1;om2];
   phi=[phi1;phi2];
   K=[K2;K1];
   Xo=[Xo1;Xo2];
   Zo=[Yo1;Yo2];
   Yo=[Zo1;Zo2];
   
   %Initial approximate values for new points 2,5,7 and 8
   Points=1000*[2.607185135191686  0.437760632224551  2.082912629762741
           3.333021059933281 0.419235848776470 2.141931670179496
           3.015317225017418  0.644350343433182 2.011771633315084
          3.013692895343422 0.435293459922647 2.177912133052440];

       Xp=Points(:,1);
       Yp=Points(:,2);
       Zp=Points(:,3);
    %%concatanation of control points with approximate values
       X=[Xc;Xp];
       Z=[Yc;Yp];
       Y=[Zc;Zp];
   
 %%
  %%%%%%FOCAL LENGTHS%%%%%%%%%%5
   f1=165.89;
   f2=165.77;
   f=[f1;f2];
 %%
%%%%%%%%%%%%%%%Customized Matrices%%%%%%%%%%%%%%
% MM=[1;1;2;2;3;3;4;4;5;5;6;6;7;7;8;8;9;9;10;10;11;11;12;12;13;13;14;14;15;15;16;16];%%pick photo coordinates%%
M=[1;1;9;9;2;2;10;10;3;3;11;11;4;4;12;12;5;5;13;13;6;6;14;14;7;7;15;15;8;8;16;16];
% GG=[1;1;5;5;2;2;3;3;6;6;4;4;7;7;8;8;1;1;5;5;2;2;3;3;6;6;4;4;7;7;8;8];%%%pick control points%%%
G=[1;1;1;1;5;5;5;5;2;2;2;2;3;3;3;3;6;6;6;6;4;4;4;4;7;7;7;7;8;8;8;8];
% HH=[1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2];%%%pick approximmate initial exterior orientation parameters Xo,Yo,Zo %%%
H=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];
% VV=[1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2];%%%pick approximate rotational elements om,phi,k %%%
V=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];
% PP=[1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2];%%%pick the focal lengths %%%%
P=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];
R=[1;1;2;2;3;3;4;4;5;5;6;6;7;7;8;8;9;9;10;10;11;11;12;12;13;13;14;14;15;15;16;16];%%%places the values in correct rows %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TT=[1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2];%%%%pick focal length%%%
T=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];
% BB=[1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2];%%%%pick w,phi,kappa%%%
B=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];
% DD=[2;2;5;5;7;7;8;8;10;10;13;13;15;15;16;16];%%%pick photo coordinates%%
D=[2;2;10;10;5;5;13;13;7;7;15;15;8;8;16;16];
EE=[1;1;2;2;3;3;4;4;1;1;2;2;3;3;4;4];%%place right column%%%
E=[1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4];
%  S=[1;1;4;4;7;7;10;10;1;1;4;4;7;7;10;10];%%place right column%%%
CC=[3;3;9;9;13;13;15;15;19;19;25;25;29;29;31;31];%%%place in correct rows%%
C=[5;5;7;7;17;17;19;19;25;25;27;27;29;29;31;31];
UU=[1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4];%%control
% S=[1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2];%%pick exterior%%
% TT=[1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2];%%%%pick focal length%%%
S=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];%%Pick exterior
 %%
 
%DIFFERENTIALS
%m is photo coordinte
%g is control point
%h is exterior approx
%n for om,phi
%p for f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%omega differentials%%
for ii=1:1;

for i=1:32;
    m=M(i);
    g=G(i);
    h=H(i);
    n=V(i);
    p=P(i);
    q=R(i);
    i=2*q-1;
    j=P(i);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 r11=cos(phi(n))*cos(K(n));
 r12=((cos(om(n))*sin(K(n)))+(sin(om(n))*sin(phi(n))*cos(K(n))));
 r13=(sin(om(n))*sin(K(n)))-(cos(om(n))*sin(phi(n))*cos(K(n)));
 r21=-cos(phi(n))*sin(K(n));
 r22=(cos(om(n))*cos(K(n)))-(sin(om(n))*sin(phi(n))*sin(K(n)));
 r23=(sin(om(n))*cos(K(n)))+(cos(om(n))*sin(phi(n))*sin(K(n)));
 r31=sin(phi(n));
 r32=-sin(om(n))*cos(phi(n));
 r33=cos(om(n))*cos(phi(n));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%
    %%%dr11,dr12,dr13 differentials %%%%
      %%%dw%%
 dr12dw=(-sin(om(n))*sin(K(n)))+(cos(om(n))*sin(phi(n))*cos(K(n)));
 dr13dw=(cos(om(n))*sin(K(n)))+(sin(om(n))*sin(phi(n))*cos(K(n)));
     %%%dphi%%
 dr11dphi=(-sin(phi(n))*cos(K(n)));
 dr12dphi=(sin(om(n))*cos(phi(n))*cos(K(n)));
 dr13dphi=(-cos(om(n))*cos(phi(n))*cos(K(n)));
     %%%dkappa%%
 dr11dK=(cos(phi(n))*sin(K(n)));
 dr12dK=(cos(om(n))*cos(K(n)))-(sin(om(n))*sin(phi(n))*sin(K(n)));
 dr13dK=(sin(om(n))*cos(K(n)))+(cos(om(n))*sin(phi(n))*sin(K(n)));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%dr21,dr22,dr23 differntials %%%%
        %%dw%%
  dr21dw=0;
  dr22dw=(-sin(om(n))*cos(K(n)))-(cos(om(n))*sin(phi(n))*sin(K(n)));
  dr23dw=(cos(om(n))*cos(K(n)))-(sin(om(n))*sin(phi(n))*sin(K(n)));
        %%dphi%%
  dr21dphi=sin(phi(n))*sin(K(n));
  dr22dphi=-sin(om(n))*cos(phi(n))*sin(K(n));
  dr23dphi=cos(om(n))*cos(phi(n))*sin(K(n));
       %%dkappa%%
  dr21dK=-cos(phi(n))*cos(K(n));
  dr22dK=(-cos(om(n))*sin(K(n)))-(sin(om(n))*sin(phi(n))*cos(K(n)));
  dr23dK=(-sin(om(n))*sin(K(n)))+(cos(om(n))*sin(phi(n))*cos(K(n)));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%dr31,dr32,dr33 differentials %%%%%%
       %%dw%%
  dr32dw=-cos(om(n))*cos(phi(n));
  dr33dw=-sin(om(n))*cos(phi(n));
       %%dphi%%
  dr31dphi=cos(phi(n));
  dr32dphi=sin(om(n))*sin(phi(n));
  dr33dphi=-cos(om(n))*sin(phi(n));
 %% 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%Xbar,Ybar,Zbar%%%%%%
 Xbar=(r11*(X(g)-Xo(h))+r12*(Z(g)-Zo(h))-r13*(Y(g)-Yo(h)));
 Ybar=((r21*(X(g)-Xo(h))+r22*(Z(g)-Zo(h))-r23*(Y(g)-Yo(h))));
 Zbar=((r31*(X(g)-Xo(h))+r32*(Z(g)-Zo(h))-r33*(Y(g)-Yo(h))));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%
 
 %%omega differentials%%
A(i,j)=(-f(p)*((((dr12dw*(Z(g)-Zo(h)))-((dr13dw*(Y(g)-Yo(h)))))*Zbar)-(((dr32dw*(Z(g)-Zo(h)))-(dr33dw*(Y(g)-Yo(h))))*Xbar)))/Zbar^2;
A(i+1,j)=(-f(p)*((((dr21dw*(X(g)-Xo(h)))+(dr22dw*(Z(g)-Zo(h)))-((dr23dw*(Y(g)-Yo(h)))))*Zbar)-(((dr32dw*(Z(g)-Zo(h)))-(dr33dw*(Y(g)-Yo(h))))*Ybar)))/Zbar^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%phi differentials%%
A(i,2+j)=(-f(p)*((((dr11dphi*(X(g)-Xo(h)))+(dr12dphi*(Z(g)-Zo(h)))-(dr13dphi*(Y(g)-Yo(h))))*Zbar)-(((dr31dphi*(X(g)-Xo(h)))+(dr32dphi*(Z(g)-Zo(h)))-(dr33dphi*(Y(g)-Yo(h))))*Xbar)))/Zbar^2;
A(i+1,2+j)=(-f(p)*((((dr22dphi*(Z(g)-Zo(h)))-(dr23dphi*(Y(g)-Yo(h))))*Zbar)-(((dr31dphi*(X(g)-Xo(h)))+(dr32dphi*(Z(g)-Zo(h)))-(dr33dphi*(Y(g)-Yo(h))))*Ybar)))/Zbar^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%kappa differentials%%
A(i,4+j)=(-f(p)*(((dr11dK*(X(g)-Xo(h)))+(dr12dK*(Z(g)-Zo(h)))-(dr13dK*(Y(g)-Yo(h))))*Zbar))/Zbar^2;
A(i+1,4+j)=(-f(p)*(((dr21dK*(X(g)-Xo(h)))+(dr22dK*(Z(g)-Zo(h)))-(dr23dK*(Y(g)-Yo(h))))*Zbar))/Zbar^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Xo differentials %%%
A(i,6+j)=(-f(p)*(-r11*Zbar+r31*Xbar))/Zbar^2;
A(i+1,6+j)=(-f(p)*(-r21*Zbar+r31*Ybar))/Zbar^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Yo differentials%%%
A(i,8+j)=(-f(p)*(r13*Zbar-r33*Xbar))/Zbar^2;
A(i+1,8+j)=(-f(p)*(r23*Zbar-r33*Ybar))/Zbar^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Zo differentials%%
A(i,10+j)=(-f(p)*(-r12*Zbar+r32*Xbar))/Zbar^2;
A(i+1,10+j)=(-f(p)*(-r22*Zbar+r32*Ybar))/Zbar^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%A matrix of coordinates of the points to be obtained%%
for i=1:16;
u=UU(i); %%%%pick control point %%%
s=S(i);%%%%%%pick exterior orientation %%%    
t=T(i); %%%%focal length%%
b=B(i);%%%%pick w,phi,kappa%%%
d=D(i);
e=E(i);
j2=3*e;
c=C(i);
i=c;


 r11=cos(phi(b))*cos(K(b));
 r12=((cos(om(b))*sin(K(b)))+(sin(om(b))*sin(phi(b))*cos(K(b))));
 r13=(sin(om(b))*sin(K(b)))-(cos(om(b))*sin(phi(b))*cos(K(b)));
 r21=-cos(phi(b))*sin(K(b));
 r22=(cos(om(b))*cos(K(b)))-(sin(om(b))*sin(phi(b))*sin(K(b)));
 r23=(sin(om(b))*cos(K(b)))+(cos(om(b))*sin(phi(b))*sin(K(b)));
 r31=sin(phi(b));
 r32=-sin(om(b))*cos(phi(b));
 r33=cos(om(b))*cos(phi(b));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 dr12dw=(-sin(om(b))*sin(K(b)))+(cos(om(b))*sin(phi(b))*cos(K(b)));
 dr13dw=(cos(om(b))*sin(K(b)))+(sin(om(b))*sin(phi(b))*cos(K(b)));
     %%%dphi%%
 dr11dphi=(-sin(phi(b))*cos(K(b)));
 dr12dphi=(sin(om(b))*cos(phi(b))*cos(K(b)));
 dr13dphi=(-cos(om(b))*cos(phi(b))*cos(K(b)));
     %%%dkappa%%
 dr11dK=(cos(phi(b))*sin(K(b)));
 dr12dK=(cos(om(b))*cos(K(b)))-(sin(om(b))*sin(phi(b))*sin(K(b)));
 dr13dK=(sin(om(b))*cos(K(b)))+(cos(om(b))*sin(phi(b))*sin(K(b)));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%dr21,dr22,dr23 differntials %%%%
        %%dw%%
  dr21dw=0;
  dr22dw=(-sin(om(b))*cos(K(b)))-(cos(om(b))*sin(phi(b))*sin(K(b)));
  dr23dw=(cos(om(b))*cos(K(b)))-(sin(om(b))*sin(phi(b))*sin(K(b)));
        %%dphi%%
  dr21dphi=sin(phi(b))*sin(K(b));
  dr22dphi=-sin(om(b))*cos(phi(b))*sin(K(b));
  dr23dphi=cos(om(b))*cos(phi(b))*sin(K(b));
       %%dkappa%%
  dr21dK=-cos(phi(b))*cos(K(b));
  dr22dK=(-cos(om(b))*sin(K(b)))-(sin(om(b))*sin(phi(b))*cos(K(b)));
  dr23dK=(-sin(om(b))*sin(K(b)))+(cos(om(b))*sin(phi(b))*cos(K(b)));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%dr31,dr32,dr33 differentials %%%%%%
       %%dw%%
  dr32dw=-cos(om(b))*cos(phi(b));
  dr33dw=-sin(om(b))*cos(phi(b));
       %%dphi%%
  dr31dphi=cos(phi(b));
  dr32dphi=sin(om(b))*sin(phi(b));
  dr33dphi=-cos(om(b))*sin(phi(b));

             %%%%%%Xbar,Ybar,Zbar%%%%%%
 Xbar=(r11*(X(u)-Xo(s))+r12*(Z(u)-Zo(s))-r13*(Y(u)-Yo(s)));
 Ybar=((r21*(X(u)-Xo(s))+r22*(Z(u)-Zo(s))-r23*(Y(u)-Yo(s))));
 Zbar=((r31*(X(u)-Xo(s))+r32*(Z(u)-Zo(s))-r33*(Y(u)-Yo(s))));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(i,10+j2)=(-f(t)*(r11*Zbar-r31*Xbar))/Zbar^2;
A(i+1,10+j2)=(-f(t)*(r21*Zbar-r31*Ybar))/Zbar^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Ya differentials%%
A(i,11+j2)=(-f(t)*(-r13*Zbar+r33*Xbar))/Zbar^2;
A(i+1,11+j2)=(-f(t)*(-r23*Zbar+r33*Ybar))/Zbar^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Za differentials %%
A(i,12+j2)=(-f(t)*(r12*Zbar-r32*Xbar))/Zbar^2;
A(i+1,12+j2)=(-f(t)*(r22*Zbar-r32*Ybar))/Zbar^2;

end
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Matrix of constants %%%%%%%%%%%%%%%%%%%%
for i=1:32;
     m=M(i);
    g=G(i);
    h=H(i);
    n=V(i);
    p=P(i);
    q=R(i);
    i=2*q-1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 r11=cos(phi(n))*cos(K(n));
 r12=((cos(om(n))*sin(K(n)))+(sin(om(n))*sin(phi(n))*cos(K(n))));
 r13=(sin(om(n))*sin(K(n)))-(cos(om(n))*sin(phi(n))*cos(K(n)));
 r21=-cos(phi(n))*sin(K(n));
 r22=(cos(om(n))*cos(K(n)))-(sin(om(n))*sin(phi(n))*sin(K(n)));
 r23=(sin(om(n))*cos(K(n)))+(cos(om(n))*sin(phi(n))*sin(K(n)));
 r31=sin(phi(n));
 r32=-sin(om(n))*cos(phi(n));
 r33=cos(om(n))*cos(phi(n));
 
    
 %% 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %%%%%%Xbar,Ybar,Zbar%%%%%%
 Xbar=(r11*(X(g)-Xo(h))+r12*(Z(g)-Zo(h))-r13*(Y(g)-Yo(h)));
 Ybar=((r21*(X(g)-Xo(h))+r22*(Z(g)-Zo(h))-r23*(Y(g)-Yo(h))));
 Zbar=((r31*(X(g)-Xo(h))+r32*(Z(g)-Zo(h))-r33*(Y(g)-Yo(h))));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lo(i,1)=(-f(p)*Xbar)/Zbar;
Lo(i+1,1)=(-f(p)*Ybar)/Zbar;

end

%%
dL=L-Lo;
N=A'*A; %%%%%Normal equation matrix %%%%%%%
Qxx=inv(N);%%%%%%%% coffactor matrix %%%%%%
d=A'*dL;   %%%%% absolute vector %%%%%%%%
delta=Qxx*d; %%%%% corrections %%%%%%%%

%%
%%%%updating initial values %%%
om1=om1+delta(1,1);
om2=om2+delta(2,1);
phi1=phi1+delta(3,1);
phi2=phi2+delta(4,1);
K1=K1+delta(5,1);
K2=K2+delta(6,1);
Xo1=Xo1+delta(7,1);
Xo2=Xo2+delta(8,1);
Yo1=Yo1+delta(9,1);
Yo2=Yo2+delta(10,1);
Zo1=Zo1+delta(11,1);
Zo2=Zo2+delta(12,1);
Xp(1,1)=Xp(1,1)+delta(13,1);
Yp(1,1)=Yp(1,1)+delta(14,1);
Zp(1,1)=Zp(1,1)+delta(15,1);
Xp(2,1)=Xp(2,1)+delta(16,1);
Yp(2,1)=Yp(2,1)+delta(17,1);
Zp(2,1)=Zp(2,1)+delta(18,1);
Xp(3,1)=Xp(3,1)+delta(19,1);
Yp(3,1)=Yp(3,1)+delta(20,1);
Zp(3,1)=Zp(3,1)+delta(21,1);
Xp(4,1)=Xp(4,1)+delta(22,1);
Yp(4,1)=Yp(4,1)+delta(23,1);
Zp(4,1)=Zp(4,1)+delta(24,1);
 
%%%%%%%%%%%solution%%%%%%%%%%%%%%%%%%%
end
 SOL=[om1;om2;phi1;phi2;K1;K2;Xo1;Xo2;Yo1;Yo2;Zo1;Zo2;Xp(1,1);Yp(1,1);Zp(1,1);Xp(2,1);Yp(2,1);Zp(2,1);Xp(3,1);Yp(3,1);Zp(3,1);Xp(4,1);Yp(4,1);Zp(4,1)];
