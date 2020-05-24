%@ Author
 clear;
 format long
 %Data given to compute parameters
 %Photo coordinates
 Photo=[-64.097	72.924
        64.017	77.183
        -23.087	-39.191
        64.353	-37.129];
    %Control points coordinates
  Control=[2294.96	4735.29	
	       2274.05	5263.62	
	       2736.69	4928.94	
	       2734.89	5269.75];
     %Elevations of the Ground points     
  Z=[1236.46;1242.56;1313.25;1274.9];
   
  havg=1239.51;%average Elevations of A and B
%   Extracting coordinates from the given data
  x=Photo(:,1);
  y=Photo(:,2);
  X=Control(:,1);
  Y=Control(:,2);
%   TC is a column matrix that will enable the loops to form a matrix.
  TC= [1; 1; 2; 2; 3; 3; 4; 4]; 
  f=0.11205;%camera focal length.
%   Computation of the initial conditions.
%Computing Kappa
 W=atan((y(2,1)-y(1,1))/(x(2,1)-x(1,1)));
 U=atan((Y(2,1)-Y(1,1))/(X(2,1)-X(1,1)));
 K=W-U;%kappa
 %Computing magnitudes of AB and ab to obtain scale.
 T=sqrt(((y(2,1)-y(1,1))^2)+((x(2,1)-x(1,1))^2));
 S=sqrt(((Y(2,1)-Y(1,1))^2)+((X(2,1)-X(1,1))^2));
 Scale=(T/1000)/S;
 %furtther initial values.
 Xo=X(1,1)-((((x(1,1)*cos(-K))+(y(1,1)*sin(-K)))/Scale)/1000);
 Yo=Y(1,1)-((((-x(1,1)*sin(-K))+(y(1,1)*cos(-K)))/Scale)/1000);
 Zo=((f+havg*Scale)/Scale);
 om=0;
 phi=0;
 %Formation of the Rotation matrix.
 r11=cos(phi)*cos(K);
 r12=((cos(om)*sin(K))+(sin(om)*sin(phi)*cos(K)));
 r13=(sin(om)*sin(K))-(cos(om)*sin(phi)*cos(K));
 r21=-cos(phi)*sin(K);
 r22=(cos(om)*cos(K))-(sin(om)*sin(phi)*sin(K));
 r23=(sin(om)*cos(K))+(cos(om)*sin(phi)*sin(K));
 r31=sin(phi);
 r32=-sin(om)*cos(phi);
 r33=cos(om)*cos(phi);
 
 
%Iterations loop
%Formation of A matrix.
%Computing coefficients of dXo
for i=1:5
for i=1:8;
      g=TC(i);
      i=2*g-1;
      dXo(i,1)=-x(g)*r31-f*r11;
      dXo(i+1,1)=-y(g)*r31-f*r21;
    
 end
%Computing coefficients of dYo
  for i=1:8;
    h=TC(i);
    i=2*h-1;
   dYo(i,1)=-x(h)*r32-f*r12;
   dYo(i+1,1)=-y(h)*r32-f*r22;
  end

%Computing coefficients of dZo
   for i=1:8;
    k=TC(i);
    i=2*k-1;
   dZo(i,1)=-x(k)*r33-f*r13;
   dZo(i+1,1)=-y(k)*r33-f*r23;
   end;
%Computing coefficients of domega
% parts forming the linearized equations are computed independently due to
% the length of the equation
   for i=1:8;
    m=TC(i);
    i=2*m-1;
  FX1=(x(m)*((-(cos(om)*cos(phi)*(Y(m)-Yo)))-(sin(om)*cos(phi)*(Z(m)-Zo)))/1000);
  FX2=f*(-(sin(om)*sin(K))+cos(om)*sin(phi)*cos(K)*(Y(m)-Yo)+((cos(om)*sin(K)+sin(om)*sin(phi)*cos(K))*(Z(m)-Zo)));
  dom(i,1)=FX1+FX2;
  FY1=(y(m)*((-(cos(om)*cos(phi)*(Y(m)-Yo)))-(sin(om)*cos(phi)*(Z(m)-Zo)))/1000);
  FY2=f*((-sin(om)*cos(K))-cos(om)*sin(phi)*sin(K))*(Y(m)-Yo)+(cos(om)*cos(K) -sin(om)*sin(phi)*sin(K)*(Z(m)-Zo));
 dom(i+1,1)=FY1+FY2;

   end
   %Computing coefficients of dphi
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
   %Computing coefficients of dkappa
   for i=1:8;
    p=TC(i);
    i=2*p-1;
   dkappa(i,1)=f*((-cos(phi)*sin(K)*(X(p)-Xo))+(cos(om)*cos(K)-sin(om)*sin(phi)*sin(K)*(Y(p)-Yo))+(sin(om)*cos(K)+cos(om)*sin(phi)*sin(K)*(Z(p)-Zo)));
   dkappa(i+1,1)=-f*(((cos(phi)*cos(K))*(X(p)-Xo))+((cos(om)*sin(K)+sin(om)*sin(phi)*cos(K))*(Y(p)-Yo))+(sin(om)*sin(K)-cos(om)*sin(phi)*cos(K))*(Z(p)-Zo));
  end

A=[dXo dYo dZo dom dphi dkappa];
%Computing matrix of constants Fo
   for i=1:8;
    q=TC(i);
    i=2*q-1;
    Fo(i,1)=x(q)*((r31*(X(q)-Xo))+(r32*(Y(q)-Yo))+(r33*(Z(q)-Zo)))/1000+f*((r11*(X(q)-Xo))+(r12*(Y(q)-Yo))+(r13*(Z(q)-Zo)));
    Fo(i+1,1)=y(q)*((r31*(X(q)-Xo))+(r32*(Y(q)-Yo))+(r33*(Z(q)-Zo)))/1000+f*((r21*(X(q)-Xo))+(r22*(Y(q)-Yo))+(r23*(Z(q)-Zo))) ; 
   end

N=A'*A;%Normal equation matrix
Qxx=inv(N);%Coffactor matrix
d=A'*-Fo;%Absolute vactor
xx=Qxx*d;%corrections to be effected
%updating of initial values
Xo=Xo+xx(1,1);
Yo=Yo+xx(2,1);
Zo=Zo+xx(3,1);
om=om+xx(4,1);
phi=phi+xx(5,1);
K=K+xx(6,1);
end
%Final exterior orientation parameters.
unqsoln=[Xo;Yo;Zo;om;phi;K]
% % computing least squares solution.
% N=A'*A;
% Qxx=inv(N);
% d=A'*-Fo;
% delta=Qxx*d;
% 
%  Xo=Xo+delta(1,1);
%  Yo=Yo+delta(2,1);
%  Zo=Zo+delta(3,1);
%  om=om+delta(4,1);
%  phi=phi+delta(5,1);
%  K=K+delta(6,1);
% Lsoln=[Xo;Yo;Zo;om;phi;K]








