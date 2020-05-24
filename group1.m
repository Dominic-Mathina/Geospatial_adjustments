% CODE DONE BY CHERONO ANITAH
% REG NO: F19/1705/2013
%CODE FOR ERROR OR OUTLIER CORRECTION INTRODUCED TO AN OBSERVATION
SP = [1; 1 ;1 ;1; 1; 1; 1; 1; 1; 1; 1; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 3; 3; 3; 3; 3; 3; 3; 3; 3; 3; 3; 4; 4; 4; 4; 4; 4; 4; 4; 4; 4; 4; 6; 6; 6; 6; 6;];
TP = [2;11;12;3;7;4;5;6;8;9;10;1;8;9;10;11;12;7;3;4;5;6;1;8;2;9;10;11;12;7;4;5;6;1;8;2;9;3;10;11;7;12;5;6;1;6;4;3;2];
D = [0;3.92813;26.65768;28.14169;29.46106;42.51533;79.96743;125.68087;272.63544;350.12485;379.70302;0;28.95076;114.46143;170.86868;205.16361;233.14019;236.90870;239.06410;252.26800;306.12031;365.52605;0;7.24369;10.92332;35.65047;71.45038;119.97324;195.22762;204.68251;243.91905;347.23966;383.60808;0;3.71152;9.75390;24.10203;29.54640;45.87412;66.81895;86.84505;97.16398;370.06912;388.29153;0;380.49429;132.61652;95.41328;26.15275];
STW = [1;2;3;4;5;6;7;8;9;10;11;12];
sp = [1;1;1;2;2;3];
tp = [2;3;4;3;4;4];
X = [100.1030;111.6010;122.1810;116.6920;87.6610;86.8550;129.5510;102.4480;126.6760;143.9770;145.6870;133.6100];
Y = [100.0110;109.0030;144.0130;168.0140;134.1990;106.2100;161.8670;90.1670;96.8140;115.7720;140.4290;163.0790];
S=[14.5967;49.2302;69.9973;36.5734;59.2302;24.6209];
% error of 0.008 introduced to observation 1
% dist = [(S(1,1)+0.008);S(2:6,4)]
% matrix for distances
for   i= 1:6;
        s=sp(i);
        t=tp(i);
        deltax= X(t)-X(s);
        deltay= Y(t)- Y(s);
        sc= sqrt(power (deltax,2) + power (deltay,2));
        j= 2*s-1;
        A1(i,j) = -deltax/sc;
        A1(i,j+1) = -deltay/sc;
        j= 2*t-1;
        A1(i,j) = deltax/sc;
        A1(i,j+1) = deltay/sc ; 
        ds(i) = S(i)-sc;
        dsf = ds';
end
% matrix for directions
        for i = 1:49;
       u = SP(i);
       v = TP(i);
       deltax2 = X(v)-X(u);
       deltay2 = Y(v)-Y(u);
       sc2 = sqrt(deltax2^2+deltay2^2);
       j2 = 2*u-1;
       A2(i,j2) = deltay2/(sc2^2);
       A2(i,j2+1)= -deltax2/(sc2^2);
       
       j2=2*v-1;
       A2(i,j2) = -deltay2/(sc2^2);
       A2(i,j2+1)= deltax2/(sc2^2);
        end
%         A = []
% conversion of grads to radians
       Drad = degtorad(D*360/400);
       % L matrix
       Dir=[0.006004725
0.006004761
0.006005482
0.006002594
0.006004464
0.006008522
0.00600353
-0.060042837
0.006004595
0.0060013
0.006002864
-0.000885404
-0.000895906
-0.000886103
-0.000885849
-0.000893817
-0.000890511
-0.00088666
-0.000891942
-0.000899945
-0.000890824
0.008906963
-0.002637813
-0.002628415
-0.002627926
-0.002636812
-0.002643222
-0.002639072
-0.002638145
-0.002646034
-0.00263777
-0.002636193
0.026371401
-0.002419699
-0.002418774
-0.002418716
-0.002414115
-0.002410034
-0.002415644
-0.002415444
-0.002410641
-0.002413796
-0.002422643
0.024159505
-0.014280394
0.05714302
-0.014286201
-0.014289764
-0.014286661
];
 DL = [dsf;Dir];       
% matrix for nuisance orientation parameters
 for i =1:11; j=1;
       B(i,j) =1;
       i= 12:22; j=2;
       B(i,j) =1;
       i = 23:33; j=3;
       B(i,j) =1;
       i= 34:44;  j=4;
       B(i,j) =1;
       i = 45:49; j=5;
       B(i,j) =1;
 end
 %combined A matrix
 A = [A1 zeros(6,21);A2,B];
 %weight matrix for distances
 for i = 1:6;
     j = 1:6;
 w1 = 1/0.0003^2*eye(6);
 end
 %weight matrix for directions
 for i=1:49;
     j=1:49;
 w2 = 1/((3*10^-3)*(pi*200))^2*eye(49);
 
 end
 %combined weight matrix
 W = [w2 zeros(49,6);zeros(6,49) w1];
 %the normal equation matrix
 N = A'*W*A;
 %extracting normal equation matrix for parameters and knowns
 Nxx = N(1:24,1:24);
 Nxt = N(1:24,25:29);
 Ntx = N(25:29,1:24);
 Ntt = N(25:29,25:29);
 %reduced normal equation matrix
 Nxx1 = Nxx-Nxt*inv(Ntt)*Ntx;
 n = A'*W*DL;
 nx = n(1:24);
 nt = n(25:29);
 %reduced absolute vector
 nx_ = nx-Nxt*inv(Ntt)*nt;
%  for loop of G matrix for free network
 for i=1;j=1:2:24;
      G(i,j)=1;
      i=2;j=2:2:24;
      G(i,j)=1;
       i=3;j=1:2:24;
       G(i,j)= Y;
      j=2:2:24;
      G(i,j)= -X;
 end
 %constraint matrix R
 R = ([[Nxx1;G] [G'; zeros(3,3)]]);
%  inverse of the constraint matrix
 Q = inv(R);
 T = [nx_ ;zeros(3,1)];
 Z = Q * T;
 dx = Z(1:24,1);
 dp = Z(25:27,1);
 qx = dx(1:2:24);
 qy = dx(2:2:24);
%  final X and Y coordinates after adjustments
 Xf = X +qx;
 Yf = Y + qy;
%  residual vextor
 v = DL-A*n;
%  aposteriori variance
 poste = (v'*W*v)/(55-24);
 Qxx = R(1:24,1:24);
 Exx = abs(Qxx*poste);
%  extracting respective covariance matrices
%  Ex = Exx(1:12,13:24);
%  Ey = Exx(13:24,1:12);
% %  extraction of variances from covariance matrix
%  varix = diag(Ex);
%  variy = diag(Ey);
% %  standard deviations
%  sdx = sqrt(varix);
%  sdy = sqrt(variy);
%  covx = diag(fliplr(Ex));
%  covy = diag(fliplr(Ey));
%  elements of error ellipse (a & b)
% direction of major axis for distances
%   Extracting covariance matrices associated with the 12 points
  E1 = diag(Exx);
  E2 = diag(Exx,1);
%%   
%  a for loop extracting the variances from the above matrix
 for j=1;
     i=1:2:24;
     sigx2=(E1(i));
     i1=1:2:23;
     Covarx = (E2(i1));
      i1=1:2:23;
     i=2:2:24;
     Covary = (E2(i1));
     sigy2= (E1(i));
 end
%%  
% obtaining elements of the error ellipses
 sigx=sqrt(sigx2);
 sigy=sqrt(sigy2);
 a=sqrt((0.5*(sigx2+sigy2))+sqrt(0.25*(sigx2-sigy2)).^2+Covarx.^2);
 b=sqrt((0.5*(sigx2+sigy2))-sqrt(0.25*(sigx2-sigy2)).^2+Covarx.^2);

%  direction of major axis for directions
     for i = 1:12;
         j=1;
        tanthita(i) = (2*Covarx(i))/(sigx2(i)-sigy2(i));
     end
        thita = (atan(tanthita)/2)';
             for i = 1:12;
            x=X(i);
            y=Y(i);
            grid on;
            hold on;
            plot(x,y);
     end
    xlabel('X coordinates');
    ylabel('Y coordinates');
    title('ERROR ELLIPSES');
 for i = 1:49;
     p = X(SP(i));
     g = X(TP(i));
     q = Y(SP(i));
     r = Y(TP(i));
     Xo = [p g];
     Yo = [q r];
     plot(Yo,Xo);
 end
%  contol size of ellipse
    for i = 1:12;
        p = 1:1:12;
        t = -2*pi:pi/100:2*pi;
        x = X(i)+((a(i)*sin(t))*25000);
        y = Y(i)+((b(i)*cos(t))*25000);
        plot(y,x);
        text(Y(i),X(i),num2str(i));
    end
 %plotting the error ellipse
 %activating grid   

%global model test
 GMT = v'*W*v/poste;
if GMT<48&GMT>16.1680;
    disp ('accept null hypothesis');
 else
   disp ('reject null hypothesis');
end
% computing the standard deviation of residuals
vmean = (sum(v))/55;
vsq = power(v-vmean,2);
sd = sqrt((sum(vsq))/55);
% detection of outliers by method of outliers
for jj  = 1:length(v);
    if (v(jj)<=3*sd);
        p(jj)=1;
    elseif (v(jj)>3*sd);
        p(jj)=exp(-v(jj));
    end
end
        
        

     




 
 
 

     
 
 
               
 
 
 
 

 
 