%load stand points and target points for directions
load('C:\Users\Tarquin\Documents\Book1.txt')
SP = Book1(1:49,1);
TP = Book1 (1:49,2);
Directions =Book1 (1:49,3);

%load approximate coordinates
load('C:\Users\Tarquin\Documents\coordinates.txt')
X = coordinates (1:12,1)
Y = coordinates (1:12,2)

%load observed distances
load('C:\Users\Tarquin\Documents\distances.txt')
SP1 = distances (1:6,1);
TP1= distances (1:6,2);
obsdist = distances (1:6,3);

%obtaining A matrix for distances
for i = 1:6
    s = SP1(i);
    t = TP1(i);
    dx= X(t)-X(s);
    dy= Y(t)-Y(s);
   approxdist=sqrt(dx^2+dy^2);
          j=2*s-1;
    Adist(i,j)= -dx/approxdist;
    Adist(i,j+1)= -dy/approxdist;
        j=2*t-1;
      
     Adist(i,j)= dx/approxdist;
     Adist(i,j+1)= dy/approxdist;
     Ldist = obsdist-approxdist;%Ldist is the matrix generated after calculating difference between observed and computed distances
end
%final A matrix for distances
adist=[Adist,zeros(6,21)];

   %obtaining A matrix for directions
for i = 1:49
    s = SP(i);
    t = TP(i);
     dx= X(t)-X(s);
    dy= Y(t)-Y(s);
   approxdist1=(dx^2+dy^2);
          j=2*s-1;
   Adr(i,j)= dy/approxdist1;
   Adr(i,j+1)=-dx/approxdist1;
    j=2*t-1;
    
 Adr(i,j)= -dy/approxdist1;
  Adr(i,j+1)= dx/approxdist1;
  
  if ((dx>0)&&(dy<0));
      theta(i)=atan(abs(dy/dx));
  elseif((dx<0)&&(dy>0));
      theta(i)=pi-(atan(abs(dy/dx)));
  elseif ((dx<0)&&(dy<0));
      theta(i)=atan(abs(dy/dx))+pi;
  else((dx>0)&&(dy<0));
      theta(i)=(2*pi)-atan(abs(dy/dx));
  end
  %for loop is used for the swing to be applied
  for i=1:11;j=1;
      c(i,j)=atan((Y(2,1)-Y(1,1))/(X(2,1)-X(1,1)));
        i=12:22;j=1;
      c(i,j)=atan((Y(1,1)-Y(2,1))/(X(1,1)-X(2,1)))+pi;
         i=23:33;j=1;
      c(i,j)=atan((Y(3,1)-Y(1,1))/(X(3,1)-X(1,1)))+pi;
         i=34:44;
      c(i,j)=atan((Y(4,1)-Y(1,1))/(X(4,1)-X(1,1)))+pi;
        i=45:49;
      c(i,j)=atan((Y(5,1)-Y(1,1))/(X(5,1)-X(1,1)));
        L=Directions+c;
  end
  %orienting the observations
  for ii=length(L)
      if(L(ii)>0&L(ii)<=pi)
          Or(ii)=L(ii);
          
      elseif (L(ii)<=1.5*pi&L(ii)>pi)
      Or(ii)=L(ii);
      
      elseif(L(ii)<=2*pi & L(ii)>1.5*pi)
          Or(ii)=L(ii);
      
      elseif (L(ii)>2*pi)
          Or(ii)=L(ii)-2*pi;
          
      else(L(ii)<=0)
          Or(ii)=L(ii)+2*pi;
      end
  end
  
  Theta=theta';
  Oriented=Or';
          
  %Ldr is the Lmatrix for directions
  Ldr=Oriented-Theta;
  
  for i=1:11;j=1;
      B(i,j)=1;
      i=12:22;j=2;
      B(i,j)=1;
      i=23:33;j=3;
      B(i,j)=1;
      i=34:44;j=4;
      B(i,j)=1;
      i=45:49;j=5;
      B(i,j)=1;
  end
  
end

%orientation parameter matrix
for i=1:49
    s=SP(i);
    om(i,s)=1;
end

%combining direction matrix with orientation parameter
adr=[Adr,om];

%overall A matrix
A=[adist;adr];

%weight matrix for distances
   for i=1:6;
       j=1:6;
       w=1/0.003^2*eye(6);
   end
      W=diag(w);
      
%weightmatrix for directions
   wdr=1/((3*10^-3)*(pi*200))^2*eye(49)
   
%combined weight matrix
w1=[w,zeros(6,49);zeros(49,6),wdr];

%normal matrix
N=A'*w1*A;

%reduced normal matrix
Nx=N(1:24,1:24);
Nxt=N(1:24,25:29);
Ntx=N(25:29,1:24);
Nt=N(25:29,25:29);
Nxx=Nx-(Nxt*inv(Nt)*Ntx); %Nxx is the reduced normal matrix

%combined L matrix
Lf=[Ldist;Ldr'];

%absolute vector n
n=A'*w1*Lf;

%reduced absolute vector nxx
nx=n(1:24,1);
nt=n(25:29,1);
nxx=nx-(Nxt*inv(Nt)*nt);

%G matrix
G=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0;
    0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1;
    100.011 -100.103 109.003 -111.601 144.013 -122.181 168.014 -116.692 134.199 -87.661 106.21 -86.855 161.867 -129.551 90.167 -102.448 96.814 -126.676 115.772 -143.977 140.429 -145.687 163.079 -133.61];
    
%combined square matrix of Nx and G
G1=[Nx,G';G,zeros(3)];

%correction matrix
C=G1*[nxx;zeros(3,1)];
c1=C (1:24,1);
c2=C(25:27,1);

c11=c1(1:2:24);
c12=c1(2:2:25);

%adjusted coordinates
Xad=X+c11;
Yad=Y+c12;

%V matrix
V=Lf-A*n;

%aposteriori variance
Ap=V'*w1*V;
Qxx=G1(1:24,1:24);
%covariance matrix
Exx=abs(Qxx*Ap);
Ex=Exx(1:12,13:24);
Ey=Exx(13:24,1:12);
VarX=diag(Ex);
CovarX=diag(fliplr(Ex));
VarY=diag(Ey);
CovarY=diag(fliplr(Ey));

%Standard deviation
SDX=sqrt(VarX);
SDY=sqrt(VarY);

%elements of error ellipses
a=sqrt((0.5*(VarX+VarY)+sqrt(0.25*(VarX-VarY)).^2+CovarX.^2));
b=sqrt((0.5*(VarX+VarY)-sqrt(0.25*(VarX-VarY)).^2+CovarX.^2));

%direction of major axis
for i=1:12
    j=1;
    g=CovarX(i);
    h=VarX(i);
    k=VarY(i);
    eta=(atan(2*g/(h-k))/2);
end

%plot of error ellipses
x=Xad;
y=Yad;

%activating grid
for i =12;
     x = X(i);
     y = Y(i);
    grid on 
    hold on
    plot (x,y);
end

xlabel('X Coordinates')
ylabel('Y Coordinates')
title('MUTIE F.N')

for i=1:49;
    p=X(SP(i));
    f=X(TP(i));
    q=Y(SP(i));
    r=Y(TP(i));
    x0=[p f];
    y0=[q r];
    plot (y0,x0);
   
end
for i=1:12;
    p=1:1:12;
    tt=-2*pi:pi/100:2*pi;
    x=X(i)+((a(i)*sin(tt)*50));
    y=Y(i)+((b(i)*cos(tt)*50));
    plot(y,x);
    text(Y(i),X(i),num2str(i));
end
 
 
%Global Model Test
GM=V'*w1*V
Vmean=(sum(V))/55;
Vs=power(V-Vmean,2);
Sd=sqrt((sum(Vs))/55);

for jj=1:length(Vs)
    if(Vs(jj)<=3*Sd)
        p(jj)=1;
        
    elseif (Vs(jj)>3*Sd)
        p(jj)=exp(-Vs(jj));
    end
end













        
     
    