%Obsv=xlsread('Data.xlsx');
% rdQ=Obsv(:,2);
 %rdF=Obsv(:,3);
 clear;
 format long
 Photo=[-64.097	72.924
        64.017	77.183
        -23.087	-39.191
        64.353	-37.129];
  Control=[2294.96	4735.29	
	       2274.05	5263.62	
	       2736.69	4928.94	
	       2734.89	5269.75];
           
   Z=[1236.46;1242.56;1313.25;1274.9];
%   Xo= 2567.856907155749;
%   Yo= 5015.889111413525;
%   Zo=1701.718123048403;
%   om=-0.011546195418;
%  phi=-0.057043359032;
%   K=1.593433644139;
   
  havg=1239.51;
  x=Photo(:,1);
  y=Photo(:,2);
  
  X=Control(:,1);
  Y=Control(:,2);
  TC= [1; 1; 2; 2; 3; 3; 4; 4]; 
  xa=-64.097;ya=72.924;	XA=2294.96;	YA=4735.29;	Za=1236.46;
xb=64.017;yb=77.183;	XB=2274.05;	YB=5263.62;	
f=0.11205;
 W=atan((yb-ya)/(xb-xa));
 U=atan((YB-YA)/(XB-XA));
 T=sqrt(((yb-ya)^2)+((xb-xa)^2));
 S=sqrt(((YB-YA)^2)+((XB-XA)^2));
  
  K=W-U;
  %KI=Ka*180;
 % K=KI/pi;
 Scale=(T/1000)/S;
 Xo=XA-((((xa*cos(-K))+(ya*sin(-K)))/Scale)/1000);
  Yo=YA-((((-xa*sin(-K))+(ya*cos(-K)))/Scale)/1000);
  Zo=((f+havg*Scale)/Scale);
%  Ka=89.6387;
%  K=degtorad(Kp);
  om=0;
  phi=0;
 r11=cos(phi)*cos(K);
 r12=((cos(om)*sin(K))+(sin(om)*sin(phi)*cos(K)));
 r13=(sin(om)*sin(K))-(cos(om)*sin(phi)*cos(K));
 r21=-cos(phi)*sin(K);
 r22=(cos(om)*cos(K))-(sin(om)*sin(phi)*sin(K));
 r23=(sin(om)*cos(K))+(cos(om)*sin(phi)*sin(K));
 r31=sin(phi);
 r32=-sin(om)*cos(phi);
 r33=cos(om)*cos(phi);
 
 
%   for i=1:6;
%       g=TC(i);
%       i=2*g-1;
%       dX(i,1)=-x(g)*r31-f*r11;
%       dX(i+1,1)=-y(g)*r31-f*r21;
%     
%   end


  for i=1:8;
    g=TC(i);
    i=2*g-1;
   dX2(i,1)=-x(g)*r31-f*r11;
   dX2(i+1,1)=-y(g)*r31-f*r21;
  end

  for i=1:8;
    h=TC(i);
    i=2*h-1;
   dY(i,1)=-x(h)*r32-f*r12;
   dY(i+1,1)=-y(h)*r32-f*r22;
  end


   for i=1:8;
    k=TC(i);
    i=2*k-1;
   dZ(i,1)=-x(k)*r33-f*r13;
   dZ(i+1,1)=-y(k)*r33-f*r23;
   end;

   for i=1:8;
    m=TC(i);
    i=2*m-1;
%   FX1=(x(m)*((-(cos(om)*cos(phi)*(Y(m)-Yo)))-(sin(om)*cos(phi)*(Z(m)-Zo)))/1000);
  FX1=(x(m)*((-(cos(om)*cos(phi)*(Y(m)-Yo)))-(sin(om)*cos(phi)*(Z(m)-Zo)))/1000);
  FX2=f*(-(sin(om)*sin(K))+cos(om)*sin(phi)*cos(K)*(Y(m)-Yo)+((cos(om)*sin(K)+sin(om)*sin(phi)*cos(K))*(Z(m)-Zo)));
  dom(i,1)=FX1+FX2;
  FY1=(y(m)*((-(cos(om)*cos(phi)*(Y(m)-Yo)))-(sin(om)*cos(phi)*(Z(m)-Zo)))/1000);
%   D(i,1)=(x(m)*((-(cos(om)*cos(phi)*(Y(m)-Yo)))-(sin(om)*cos(phi)*(Z(m)-Zo)))/1000)+(f*(-(sin(om)*sin(K))+cos(om)*sin(phi)*cos(K)*(Y(m)-Yo)+((cos(om)*sin(K)+sin(om)*sin(phi)*cos(K))*(Z(m)-Zo))))
%   D(i+1,1)=(x(m)*((-(cos(om)*cos(phi)*(Y(m)-Yo)))-(sin(om)*cos(phi)*(Z(m)-Zo)))/1000)+(f*(-(sin(om)*sin(K))+cos(om)*sin(phi)*cos(K)*(Y(m)-Yo)+((cos(om)*sin(K)+sin(om)*sin(phi)*cos(K))*(Z(m)-Zo))))
  FY2=f*((-sin(om)*cos(K))-cos(om)*sin(phi)*sin(K))*(Y(m)-Yo)+(cos(om)*cos(K) -sin(om)*sin(phi)*sin(K)*(Z(m)-Zo));
 dom(i+1,1)=FY1+FY2;

   end
   for i=1:8;
    n=TC(i);
    i=2*n-1;
    FXP1=(x(n)*((cos(phi)*(X(n)-Xo))+sin(om)*sin(phi)*(Y(n)-Yo)-cos(om)*sin(phi)*(Z(n)-Zo)))/1000;
    FXP2=f*(-sin(phi)*cos(K)*(X(n)-Xo)+sin(om)*cos(phi)*cos(K)*(Y(n)-Yo)-cos(om)*cos(phi)*cos(K)*(Z(n)-Zo));
    dphi(i,1)=FXP1+FXP2;
    FYP1=(y(n)*((cos(phi)*(X(n)-Xo))+sin(om)*sin(phi)*(Y(n)-Yo)-cos(om)*sin(phi)*(Z(n)-Zo)))/1000;
    FYP2=f*((sin(phi)*sin(K)*(X(n)-Xo))-(sin(om)*cos(phi)*sin(K)*(Y(n)-Yo))+cos(om)*cos(phi)*sin(K)*(Z(n)-Zo));
    dphi(i+1,1)=FYP1+FYP2;
   end
   for i=1:8;
    p=TC(i);
    i=2*p-1;
   dkappa(i,1)=f*((-cos(phi)*sin(K)*(X(p)-Xo))+(cos(om)*cos(K)-sin(om)*sin(phi)*sin(K)*(Y(p)-Yo))+(sin(om)*cos(K)+cos(om)*sin(phi)*sin(K)*(Z(p)-Zo)));
   dkappa(i+1,1)=-f*(((cos(phi)*cos(K))*(X(p)-Xo))+((cos(om)*sin(K)+sin(om)*sin(phi)*cos(K))*(Y(p)-Yo))+(sin(om)*sin(K)-cos(om)*sin(phi)*cos(K))*(Z(p)-Zo));
  end

A=[dX2 dY dZ dom dphi dkappa];

   for i=1:8;
    q=TC(i);
    i=2*q-1;
    Fo(i,1)=x(q)*((r31*(X(q)-Xo))+(r32*(Y(q)-Yo))+(r33*(Z(q)-Zo)))/1000+f*((r11*(X(q)-Xo))+(r12*(Y(q)-Yo))+(r13*(Z(q)-Zo)));
    Fo(i+1,1)=y(q)*((r31*(X(q)-Xo))+(r32*(Y(q)-Yo))+(r33*(Z(q)-Zo)))/1000+f*((r21*(X(q)-Xo))+(r22*(Y(q)-Yo))+(r23*(Z(q)-Zo))) ; 
   end


% Initial=[Xo;Yo;Zo;om;phi;Ka];
% Auniq=A(1:6,1:6);
% Luniq=Fo(1:6,1);
% 
% N=Auniq'*Auniq;
N=A'*A;
Qxx=inv(N);
d=A'*-Fo;
delta=Qxx*d;

 Xo=Xo+delta(1,1);
 Yo=Yo+delta(2,1);
 Zo=Zo+delta(3,1);
 om=om+delta(4,1);
 phi=phi+delta(5,1);
 K=K+delta(6,1);
 
% N=A'*A;%Normal equation matrix
% Qxx=inv(N);%Coffactor matrix
% d=A'*-Fo2;%Absolute vactor
% xx=Qxx*d;%corrections to be effected
% %updating of initial values
% Xo=Xo+xx(1,1);
% Yo=Yo+xx(2,1);
% Zo=Zo+xx(3,1);
% om=om+xx(4,1);
% phi=phi+xx(5,1);
% K=K+xx(6,1);
% % Final solution by least squares

Lsoln=[Xo;Yo;Zo;om;phi;K]









