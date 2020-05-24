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
 x1=Photo1(:,1);
 y1=Photo1(:,2);

Photo2=[-77.425 39.310          %%Photo coordinates of photo2     
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
        3310.35 816.31 1986.42
        3263.10 245.96 2129.39];
    Xc=CTRL(:,1);
    Yc=CTRL(:,2);
    Zc=CTRL(:,3);
    



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
   %%%%%%%%%%concatanate these values%%%%%%%%%%%%
   om=[om1;om2];
   phi=[phi1;phi2];
   K=[K2;K1];
   Xo=[Xo1;Xo2];
   Yo=[Yo1;Yo2];
   Zo=[Zo1;Zo2];
   
   %Initial approximate values for new points 2,5,7 and 8
   Points=[2606 439 2080
           3320 488 2107
           3015 645 2008
           3013 424 2188];
       Xp=Points(:,1);
       Yp=Points(:,2);
       Zp=Points(:,3);
    %%concatanation of control points with approximate values
       X=[Xc;Xp];
       Y=[Yc;Yp];
       Z=[Zc;Zp];
   
 %%
  %%%%%%FOCAL LENGTHS%%%%%%%%%%5
   f1=165.89;
   f2=165.77;
   f=[f1;f2];
   
   %%%%%%%%%%%%%%%Customized Column Matrices%%%%%%%%%%%%%%
  %%%%form A matrix containing the Exterior orientation Parameters 1:12%%%%
M=[1;1;9;9;2;2;10;10;3;3;11;11;4;4;12;12;5;5;13;13;6;6;14;14;7;7;15;15;8;8;16;16];%%pick photo coordinates%%
G=[1;1;1;1;5;5;5;5;2;2;2;2;3;3;3;3;6;6;6;6;4;4;4;4;7;7;7;7;8;8;8;8];%%%pick control points%%%
H=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];%%%pick approximmate initial exterior orientation parameters Xo,Yo,Zo %%%
V=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];%%%pick approximate rotational elements om,phi,k %%%
P=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];%%%pick the focal lengths %%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:32;
    m=M(i);%%pick photo coordinates%%
    g=G(i);%%%pick control points%%%
    h=H(i);%%%pick approximmate initial exterior orientation parameters Xo,Yo,Zo %%%
    n=V(i);%%%pick approximate rotational elements om,phi,k %%%
    p=P(i);%%%pick the focal lengths %%%%
    q=R(i);
    i=2*q-1;
    j=P(i);
   %%%%%%%%%%%%%%%% 
 r11=cos(phi(n))*cos(K(n));
 r12=((cos(om(n))*sin(K(n)))+(sin(om(n))*sin(phi(n))*cos(K(n))));
 r13=(sin(om(n))*sin(K(n)))-(cos(om(n))*sin(phi(n))*cos(K(n)));
 r21=-cos(phi(n))*sin(K(n));
 r22=(cos(om(n))*cos(K(n)))-(sin(om(n))*sin(phi(n))*sin(K(n)));
 r23=(sin(om(n))*cos(K(n)))+(cos(om(n))*sin(phi(n))*sin(K(n)));
 r31=sin(phi(n));
 r32=-sin(om(n))*cos(phi(n));
 r33=cos(om(n))*cos(phi(n));
  %%%%%%%%%%%%%%%%%%%%%
        %%%dw%%
  dr11dw=0;
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
  dr31dw=0;
  dr32dw=-cos(om(n))*cos(phi(n));
  dr33dw=-sin(om(n))*cos(phi(n));
       %%dphi%%
  dr31dphi=cos(phi(n));
  dr32dphi=sin(om(n))*sin(phi(n));
  dr33dphi=-cos(om(n))*sin(phi(n));
     %%%%dK%%%
  dr31dK=0;
  dr32dK=0;
  dr33dK=0;
  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%dw%%
 A(i,j)=x(m)*((dr31dw*(X(g)-Xo(h))+dr32dw*(Y(g)-Yo(h))+dr33dw*(Z(g)-Zo(h))))+f(p)*(dr11dw*(X(g)-Xo(h))+dr12dw*(Y(g)-Yo(h))+dr13dw*(Z(g)-Zo(h)));
 A(i+1,j)=y(m)*((dr31dw*(X(g)-Xo(h))+dr32dw*(Y(g)-Yo(h))+dr33dw*(Z(g)-Zo(h))))+f(p)*(dr21dw*(X(g)-Xo(h))+dr22dw*(Y(g)-Yo(h))+dr23dw*(Z(g)-Zo(h)));
    
 %%%%%%dphi%%%
 A(i,2+j)=x(m)*((dr31dphi*(X(g)-Xo(h))+dr32dphi*(Y(g)-Yo(h))+dr33dphi*(Z(g)-Zo(h))))+f(p)*(dr11dphi*(X(g)-Xo(h))+dr12dphi*(Y(g)-Yo(h))+dr13dphi*(Z(g)-Zo(h)));
 A(i+1,2+j)=y(m)*((dr31dphi*(X(g)-Xo(h))+dr32dphi*(Y(g)-Yo(h))+dr33dphi*(Z(g)-Zo(h))))+f(p)*(dr21dphi*(X(g)-Xo(h))+dr22dphi*(Y(g)-Yo(h))+dr23dphi*(Z(g)-Zo(h)));
    
 %%%%%%dK%%%  
    
 A(i,4+j)=x(m)*((dr31dK*(X(g)-Xo(h))+dr32dK*(Y(g)-Yo(h))+dr33dK*(Z(g)-Zo(h))))+f(p)*(dr11dK*(X(g)-Xo(h))+dr12dK*(Y(g)-Yo(h))+dr13dK*(Z(g)-Zo(h)));
 A(i+1,4+j)=y(m)*((dr31dK*(X(g)-Xo(h))+dr32dK*(Y(g)-Yo(h))+dr33dK*(Z(g)-Zo(h))))+f(p)*(dr21dK*(X(g)-Xo(h))+dr22dK*(Y(g)-Yo(h))+dr23dK*(Z(g)-Zo(h))); 
 
 %%%%%%%%%dXo%%
 A(i,6+j)=-r31*x(m)-r11*f(p);
 A(i+1,6+j)=-r31*x(m)-r21*f(p);
 
 %%%%%%%dYo%%%%
 A(i,8+j)=-r32*x(m)-r12*f(p);
 A(i+1,8+j)=-r32*x(m)-r22*f(p);
 
 %%%%%%%%%dZo%%%%%%%%%
 A(i,10+j)=-r33*x(m)-r13*f(p);
 A(i+1,10+j)=-r33*x(m)-r23*f(p);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  for i=1:16  
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
 
  %%%%XA%%
 A(i,10+j2)=r31*x(m)+r11*f(p);
 A(i+1,10+j2)=r31*x(m)+r21*f(p);
 
 %%%%%%%dYA%%%%
 A(i,11+j2)=r32*x(m)+r12*f(p);
 A(i+1,11+j2)=r32*x(m)+r22*f(p);
 
 %%%%%%%%%dZA%%%%%%%%%
 A(i,12+j2)=r33*x(m)+r13*f(p);
 A(i+1,12+j2)=r33*x(m)+r23*f(p);
    
    
  end   
      
end

for i=1:32
    m=M(i);%%pick photo coordinates%%
    g=G(i);%%%pick control points%%%
    h=H(i);%%%pick approximmate initial exterior orientation parameters Xo,Yo,Zo %%%
    n=V(i);%%%pick approximate rotational elements om,phi,k %%%
    p=P(i);%%%pick the focal lengths %%%%
    q=R(i);
    i=2*q-1;
    j=P(i);
   %%%%%%%%%%%%%%%% 
 r11=cos(phi(n))*cos(K(n));
 r12=((cos(om(n))*sin(K(n)))+(sin(om(n))*sin(phi(n))*cos(K(n))));
 r13=(sin(om(n))*sin(K(n)))-(cos(om(n))*sin(phi(n))*cos(K(n)));
 r21=-cos(phi(n))*sin(K(n));
 r22=(cos(om(n))*cos(K(n)))-(sin(om(n))*sin(phi(n))*sin(K(n)));
 r23=(sin(om(n))*cos(K(n)))+(cos(om(n))*sin(phi(n))*sin(K(n)));
 r31=sin(phi(n));
 r32=-sin(om(n))*cos(phi(n));
 r33=cos(om(n))*cos(phi(n));
  %%%%%%%%%%%%%%%%%%%%%
        %%%dw%%
  dr11dw=0;
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
  dr31dw=0;
  dr32dw=-cos(om(n))*cos(phi(n));
  dr33dw=-sin(om(n))*cos(phi(n));
       %%dphi%%
  dr31dphi=cos(phi(n));
  dr32dphi=sin(om(n))*sin(phi(n));
  dr33dphi=-cos(om(n))*sin(phi(n));
     %%%%dK%%%
  dr31dK=0;
  dr32dK=0;
  dr33dK=0;%%%%%%%%%%%%%%
 %%
L(i,1)=x(m)*(r31*(X(g)-Xo(h))+r32*(Y(g)-Yo(h))+r33*(Z(g)-Zo(h)))+f(p)*(r11*(X(g)-Xo(h))+r12*(Y(g)-Yo(h))+r13*(Z(g)-Zo(h)));
L(i+1,1)=x(m)*(r31*(X(g)-Xo(h))+r32*(Y(g)-Yo(h))+r33*(Z(g)-Zo(h)))+f(p)*(r21*(X(g)-Xo(h))+r22*(Y(g)-Yo(h))+r23*(Z(g)-Zo(h)));
  %%%%%%%%%%%%%%%%%
  
    
end



N=A'*A; %%%%%Normal equation matrix %%%%%%%
Qxx=inv(N);%%%%%%%% coffactor matrix %%%%%%
d=A'*-L;   %%%%% absolute vector %%%%%%%%
delta=Qxx*d; 
%%%%% corrections %%%%%%%%
%%
%%%updating initial values %%%

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


SOLN=[om1;om2;phi1;phi2;K1;K2;Xo1;Xo2;Yo1;Yo2;Zo1;Zo2;Xp(1,1);Yp(1,1);Zp(1,1);Xp(2,1);Yp(2,1);Zp(2,1);Xp(3,1);Yp(3,1);Zp(3,1);Xp(4,1);Yp(4,1);Zp(4,1)];