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
    

for ii=1:10;
 %Initial approximate values for photo1
   Xo1=2480;
   Yo1=559;
   Zo1=378;
   om1=0;
   phi1=0;
   K1=0;
   %Initial approximate values for photo2

   Xo2=3033;
   Yo2=469;
   Zo2=377;
   om2=0;
   phi2=0;
   K2=0;
   %%%%%%%%%%concatanate these values%%%%%%%%%%%%
   om=[om1;om2];
   phi=[phi1;phi2];
   K=[K2;K1];
   Xo=[Xo1;Xo2];
   Zo=[Yo1;Yo2];
   Yo=[Zo1;Zo2];
   
   %Initial approximate values for new points 2,5,7 and 8
   Points=[2606 439 2088
           3320 488 2107
           3015 645 2008
           3013 424 2188];
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
%%%%%%%%%%%%%%%Customized Column Matrices%%%%%%%%%%%%%%
  %%%%form A matrix containing the Exterior orientation Parameters 1:12%%%%
MM=[1;1;9;9;2;2;10;10;3;3;11;11;4;4;12;12;5;5;13;13;6;6;14;14;7;7;15;15;8;8;16;16];%%pick photo coordinates%%
M=[1;1;9;9;3;3;11;11;4;4;12;12;6;6;14;14];
GG=[1;1;1;1;5;5;5;5;2;2;2;2;3;3;3;3;6;6;6;6;4;4;4;4;7;7;7;7;8;8;8;8];%%%pick control points%%%
G=[1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4];
HH=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];%%%pick approximmate initial exterior orientation parameters Xo,Yo,Zo %%%
H=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];
VV=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];%%%pick approximate rotational elements om,phi,k %%%
V=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];
PP=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];%%%pick the focal lengths %%%%
P=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];
RR=[1;1;2;2;3;3;4;4;5;5;6;6;7;7;8;8;9;9;10;10;11;11;12;12;13;13;14;14;15;15;16;16];%%%places the values in correct rows %%%
R=[1;1;2;2;3;3;4;4;5;5;6;6;7;7;8;8];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%form A matrix containing the new points;column 13:24%%%%%%%%%%%
T=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];%%%%pick focal length%%%
B=[1;1;2;2;1;1;2;2;1;1;2;2;1;1;2;2];%%%%pick w,phi,kappa%%%
D=[2;2;10;10;5;5;13;13;7;7;15;15;8;8;16;16];%%%pick photo coordinates%%
E=[1;1;1;1;4;4;4;4;7;7;7;7;10;10;10;10];%%place right column%%%
C=[5;5;7;7;17;17;19;19;25;25;27;27;29;29;31;31];%%%place in correct rows%%

 %%

for i=1:16;
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
 %%omega differentials%%
A1(i,j)=x(m)*((-cos(om(n))*cos(phi(n))*(Z(g)-Zo(h)))+(sin(om(n))*cos(phi(n)))*(Y(g)-Yo(h)))+f(p)*(((-sin(om(n))*sin(K(n)))+(cos(om(n))*sin(phi(n))*cos(K(n))))*(Z(g)-Zo(h))-((cos(om(n))*sin(K(n)))+(sin(om(n))*sin(phi(n))*cos(K(n))))*(Y(g)-Yo(h)));
A1(i+1,j)=y(m)*((-cos(om(n))*cos(phi(n))*(Z(g)-Zo(h)))+(sin(om(n))*cos(phi(n))*(Y(g)-Yo(h))))-f(p)*(((sin(om(n))*cos(K(n))+cos(om(n))*sin(phi(n))*sin(K(n)))*(Z(g)-Zo(h)))+(cos(om(n))*cos(K(n))-sin(om(n))*sin(phi(n))*sin(K(n)))*(Y(g)-Yo(h)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%phi differentials%%
A1(i,2+j)=x(m)*((cos(phi(n))*(X(g)-Xo(h)))+(sin(om(n))*sin(phi(n))*(Z(g)-Zo(h)))+(cos(om(n))*sin(phi(n))*(Y(g)-Yo(h))))+f(p)*((sin(phi(n))*cos(K(n))*(X(g)-Xo(h)))+(sin(om(n))*cos(phi(n))*cos(K(n))*(Z(g)-Zo(h)))+(cos(om(n))*cos(phi(n))*cos(K(n))*(Y(g)-Yo(h))));
A1(i+1,2+j)=y(m)*((cos(phi(n))*(X(g)-Xo(h)))+(sin(om(n))*sin(phi(n))*(Z(g)-Zo(h)))+(cos(om(n))*sin(phi(n))*(Y(g)-Yo(h))))+f(p)*((sin(phi(n))*sin(K(n))*(X(g)-Xo(h)))-(sin(om(n))*cos(phi(n))*sin(K(n))*(Z(g)-Zo(h)))-(cos(om(n))*cos(phi(n))*sin(K(n))*(Y(g)-Yo(h))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%kappa differentials%%
A1(i,4+j)=f(p)*((-cos(phi(n))*sin(K(n))*(X(g)-Xo(h)))+((cos(om(n))*cos(K(n))-sin(om(n))*sin(phi(n))*sin(K(n)))*(Z(g)-Zo(h)))-(sin(om(n))*cos(K(n))+cos(om(n))*sin(phi(n))*sin(K(n)))*(Y(g)-Yo(h)));
A1(i+1,4+j)=f(p)*((-cos(phi(n))*cos(K(n))*(X(g)-Xo(h)))+(-cos(om(n))*sin(K(n))-sin(om(n))*sin(phi(n))*cos(K(n)))*(Z(g)-Zo(h))-(-sin(om(n))*sin(K(n))+cos(om(n))*sin(phi(n))*cos(K(n))*(Y(g)-Yo(h))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Xo differentials %%%
A1(i,6+j)=-x(m)*r31-f(p)*r11;
A1(i+1,6+j)=-y(m)*r31-f(p)*r21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Yo differentials%%%
A1(i,8+j)=x(m)*r33+f(p)*r13;
A1(i+1,8+j)=y(m)*r33+f(p)*r23;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Zo differentials%%
A1(i,10+j)=-x(m)*r32-f(p)*r12;
A1(i+1,10+j)=-y(m)*r32-f(p)*r22;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%Matrix of constants %%%%%%%%%%%%%%%%%%%%
for i=1:16;
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
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L1(i,1)=x(m)*(r31*(X(g)-Xo(h))+r32*(Z(g)-Zo(h))-r33*(Y(g)-Yo(h)))+f(p)*(r11*(X(g)-Xo(h))+r12*(Z(g)-Zo(h))-r13*(Y(g)-Yo(h)));
L1(i+1,1)=y(m)*((r31*(X(g)-Xo(h))+r32*(Z(g)-Zo(h))-r33*(Y(g)-Yo(h))))+f(p)*((r21*(X(g)-Xo(h))+r22*(Z(g)-Zo(h))-r23*(Y(g)-Yo(h))));
end

N1=A1'*A1; %%%%%Normal equation matrix %%%%%%%
Qxx1=inv(N1);%%%%%%%% coffactor matrix %%%%%%
d1=A1'*-L1;   %%%%% absolute vector %%%%%%%%
delta1=Qxx1*d1; 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
om1=om1+delta1(1,1);
om2=om2+delta1(2,1);
phi1=phi1+delta1(3,1);
phi2=phi2+delta1(4,1);
K1=K1+delta1(5,1);
K2=K2+delta1(6,1);
Xo1=Xo1+delta1(7,1);
Xo2=Xo2+delta1(8,1);
Yo1=Yo1+delta1(9,1);
Yo2=Yo2+delta1(10,1);
Zo1=Zo1+delta1(11,1);
Zo2=Zo2+delta1(12,1);
%%
end
ANS=[om1;om2;phi1;phi2;K1;K2;Xo1;Xo2;Yo1;Yo2;Zo1;Zo2;];

% for i=1:16;
% t=T(i);
% b=B(i);
% d=D(i);
% e=E(i);
% j2=e;
% c=C(i);
% i=c;
%         
%  r11=cos(phi(b))*cos(K(b));
%  r12=((cos(om(b))*sin(K(b)))+(sin(om(b))*sin(phi(b))*cos(K(b))));
%  r13=(sin(om(b))*sin(K(b)))-(cos(om(b))*sin(phi(b))*cos(K(b)));
%  r21=-cos(phi(b))*sin(K(b));
%  r22=(cos(om(b))*cos(K(b)))-(sin(om(b))*sin(phi(b))*sin(K(b)));
%  r23=(sin(om(b))*cos(K(b)))+(cos(om(b))*sin(phi(b))*sin(K(b)));
%  r31=sin(phi(b));
%  r32=-sin(om(b))*cos(phi(b));
%  r33=cos(om(b))*cos(phi(b));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A2(i,0+j2)=x(d)*r31+f(t)*r11;
% A2(i+1,0+j2)=y(d)*r31+f(t)*r21;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%Ya differentials%%
% A2(i,1+j2)=-x(d)*r33-f(t)*r13;
% A2(i+1,1+j2)-y(d)*r33-f(t)*r23;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%Za differentials %%
% A2(i,2+j2)=x(d)*r32+f(t)*r12;
% A2(i+1,2+j2)=y(d)*r32+f(t)*r22;
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%Matrix of constants %%%%%%%%%%%%%%%%%%%%
% for i=1:32;
%      m=M(i);
%     g=G(i);
%     h=H(i);
%     n=V(i);
%     p=P(i);
%     q=R(i);
%     i=2*q-1;
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  r11=cos(phi(n))*cos(K(n));
%  r12=((cos(om(n))*sin(K(n)))+(sin(om(n))*sin(phi(n))*cos(K(n))));
%  r13=(sin(om(n))*sin(K(n)))-(cos(om(n))*sin(phi(n))*cos(K(n)));
%  r21=-cos(phi(n))*sin(K(n));
%  r22=(cos(om(n))*cos(K(n)))-(sin(om(n))*sin(phi(n))*sin(K(n)));
%  r23=(sin(om(n))*cos(K(n)))+(cos(om(n))*sin(phi(n))*sin(K(n)));
%  r31=sin(phi(n));
%  r32=-sin(om(n))*cos(phi(n));
%  r33=cos(om(n))*cos(phi(n));
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L(i,1)=x(m)*(r31*(X(g)-Xo(h))+r32*(Z(g)-Zo(h))-r33*(Y(g)-Yo(h)))+f(p)*(r11*(X(g)-Xo(h))+r12*(Z(g)-Zo(h))-r13*(Y(g)-Yo(h)));
% L(i+1,1)=y(m)*((r31*(X(g)-Xo(h))+r32*(Z(g)-Zo(h))-r33*(Y(g)-Yo(h))))+f(p)*((r21*(X(g)-Xo(h))+r22*(Z(g)-Zo(h))-r23*(Y(g)-Yo(h))));
% end
% %%
% A=[A1 A2];
% N=A'*A; %%%%%Normal equation matrix %%%%%%%
% Qxx=inv(N);%%%%%%%% coffactor matrix %%%%%%
% d=A'*-L;   %%%%% absolute vector %%%%%%%%
% delta=Qxx*d;
% %%%%% corrections %%%%%%%%
% %%
% %%%updating initial values %%%
% 
% om1=om1+delta(1,1);
% om2=om2+delta(2,1);
% phi1=phi1+delta(3,1);
% phi2=phi2+delta(4,1);
% K1=K1+delta(5,1);
% K2=K2+delta(6,1);
% Xo1=Xo1+delta(7,1);
% Xo2=Xo2+delta(8,1);
% Yo1=Yo1+delta(9,1);
% Yo2=Yo2+delta(10,1);
% Zo1=Zo1+delta(11,1);
% Zo2=Zo2+delta(12,1);
% Xp(1,1)=Xp(1,1)+delta(13,1);
% Yp(1,1)=Yp(1,1)+delta(14,1);
% Zp(1,1)=Zp(1,1)+delta(15,1);
% Xp(2,1)=Xp(2,1)+delta(16,1);
% Yp(2,1)=Yp(2,1)+delta(17,1);
% Zp(2,1)=Zp(2,1)+delta(18,1);
% Xp(3,1)=Xp(3,1)+delta(19,1);
% Yp(3,1)=Yp(3,1)+delta(20,1);
% Zp(3,1)=Zp(3,1)+delta(21,1);
% Xp(4,1)=Xp(4,1)+delta(22,1);
% Yp(4,1)=Yp(4,1)+delta(23,1);
% Zp(4,1)=Zp(4,1)+delta(24,1);
% 
% 
% % %%%%%%%%%%%%solution%%%%%%%%%%%%%%%%%%%
% 
%  JIBU=[om1;om2;phi1;phi2;K1;K2;Xo1;Xo2;Yo1;Yo2;Zo1;Zo2;Xp(1,1);Yp(1,1);Zp(1,1);Xp(2,1);Yp(2,1);Zp(2,1);Xp(3,1);Yp(3,1);Zp(3,1);Xp(4,1);Yp(4,1);Zp(4,1)];
% 
% % %Aa=A(1:32,1:12);
% % N = Aa'*Aa;
% % dL=Aa'*-L;
% sol =N\dLzz;