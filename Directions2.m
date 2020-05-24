
% while delX>0.02|| delX<-0.02
    
P = [1 2 0; 1 11 3.92813; 1 12 26.65768; 1 3 28.14169; 1 7 29.46106; 1 4 42.51533; 1 5 79.96743; 1 6  125.68087; ...]
 1 8 272.63544; 1 9 350.12485; 1 10 379.70302; 2 1 0; 2 8 28.95076; 2 9 114.46143; 2 10 170.86868; 2 11 205.16361; ...
 2 12 233.14019; 2 7 236.9087; 2 3 239.0641; 2 4 252.268; 2 5 306.12031; 2 6 365.52605; 3 1 0; 3 8 7.24369; 3 2 10.92332;...
 3 9 35.65047; 3 10 71.45038; 3 11 119.97324; 3 12 195.22762; 3 7 204.68251; 3 4 243.91905; 3 5 347.23966; ...
 3 6 383.60808; 4 1 0; 4 8 3.71152; 4 2 9.7539; 4 9 24.10203; 4 3 29.5464; 4 10 45.87412; 4 11 66.81895; ...
 4 7 86.84505; 4 12 97.16398; 4 5 370.06912; 4 6 388.29153; 5 1 0; 5 6 380.49429; 5 4 132.61652; 5 3 95.41328; ...
 5 2 26.1527];
 Spd = P(:,1);
 Tpd = P(:,2);
 DR = P(:,3);
STW = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12];
X = [100.103; 111.601; 122.181; 116.692; 87.661; 86.855; 129.551; 102.448; 126.676; 143.977; 145.687; 133.61];
Y = [100.011; 109.003; 144.013; 168.014; 134.199; 106.21; 161.867; 90.167; 96.814; 115.772; 140.429; 163.079];
Sp = [1; 1; 1; 2; 2; 3];
Tp = [2; 3; 4; 3; 4; 4];
DIST = [14.5967; 49.2302; 69.9973; 36.5734; 59.2302; 24.6209];
DRrad = degtorad(DR * (360/400));
% for i=1:50
% matrix for distances 
       for i = 1:6;
           s = Sp(i);
           t = Tp(i);
           dx = X(t)-X(s);
           dy = Y(t)-Y(s);
           Sc=sqrt(power(dx,2) +power(dy,2));
           j = 2*s-1 ;
           A1(i,j) = -dx/Sc;
           A1(i,j+1)= -dy/Sc;
   
           j = 2*t-1;
           A1(i,j) = dx/Sc;
           A1(i,j+1) = dy/Sc ; 
           ds(i) = DIST(i) - Sc;

       end
      dst = ds';
      
      
      for i=1:49
          k = Spd(i);
          h = Tpd(i);
          deltax = X(h)-X(k);
          deltay = Y(h)-Y(k);
          ddr(i,1)=(atan(deltay/deltax)) ;
       
      end
       ddr1=ddr(:,1);
       
% matrix for directions
   for i = 1:49;
       k = Spd(i);
       h = Tpd(i);
       deltax = X(h)-X(k);
       deltay = Y(h)-Y(k);
       ss = sqrt(deltax^2+deltay^2);
       j2 = 2*k-1;
       A2(i,j2) = deltay/(ss^2);
       A2(i,j2+1)= -deltax/(ss^2);
       
       j2=2*h-1;
       A2(i,j2) = -deltay/(ss^2);
       A2(i,j2+1)=deltax/(ss^2);
       
          j2=1;
    if ((deltax>0) &&(deltay>0));
        theta(i) = atan (abs(deltay/deltax));
    elseif ((deltax<0) && (deltay>0));
        theta(i) = pi -(atan (abs(deltay/deltax)));
    elseif ((deltax<0) && (deltay<0));
        theta(i) = atan (abs(deltay/deltax)) + pi;
    else ((deltax>0) && (deltay<0));
        theta(i) = (2*pi)-atan (abs(deltay/deltax));
    end
    %for loop for the swing to be applied
       for i=1:11;j=1;
        C(i,j) = atan((Y(2,1)-Y(1,1))/(X(2,1)-X(1,1)));
        i=12:22;j=1;
        C(i,j) = atan((Y(1,1)-Y(2,1))/(X(1,1)-X(2,1))) + pi;
        i = 23:33;j=1;
        C(i,j) = atan((Y(3,1)-Y(1,1))/(X(3,1)-X(1,1))) + pi;
        i= 34:44;
        C(i,j) = atan((Y(4,1)-Y(1,1))/(X(4,1)-X(1,1))) + pi;
        i = 45:49 ;
        C(i,j) = atan((Y(5,1)-Y(1,1))/(X(5,1)-X(1,1))) ; 
        L=DRrad+C; 
       end
         
   end 
   % An if else statetment to orient the observations

  for ii=1:length(L) 
   if (L(ii)>0 & L(ii)<=pi)
          Or(ii) = L(ii);
%          disp(L(ii))
        elseif (L(ii)<=1.5*pi & L(ii)>pi )
         Or(ii) =  L(ii);
%            disp(L(ii))
        elseif (L(ii)<=2*pi  & L(ii)>1.5*pi)
           Or(ii) = L(ii);
%          disp(L(ii))
        elseif (L(ii)>2*pi)
          Or(ii) = L(ii)-2*pi;
%           disp(L(ii)-2*pi)
        else (L(ii)<=0)
         Or(ii) = L(ii) + 2*pi;
 %           disp(L(ii) + 2*pi)
   end
 end
   
 % Computed bearings between points using provisional coordinates
   Theta1=theta';
   Dir_observed = Or';
   dL = Dir_observed - Theta1;

%difference between observed and computed bearing from excel sheet 

dLf=[dst;dL];
% The nuisance parameter for directions 
   for i=1:11;j=1;
       B(i,j) = 1;
       i=12:22;j=2;
       B(i,j) =1;
       i=23:33;j=3;
       B(i,j) =1;
       i=34:44;j=4;
       B(i,j) =1;
       i=45:49;j=5;
       B(i,j) =1;
       
   end
   
   
 for i=1;j=1:49; 
  % The weight matrix for directions
W5(i,j)=(4.25451703*10^10);
 end
 W7=diag(W5);
 
   for i=1;j=1:6; 
   % The weight matrix for distances   
W6(i,j)=(11111111.11);
  end
 W8=diag(W6);

  
 W=[W8 zeros(6,49);zeros(49,6) W7];
   At=[zeros(6,5);B];
   Ax=[A1 zeros(6,16);A2];
   Nxx=[Ax'*W*Ax];
   Nxt=[Ax'*W*At];
   Ntx=[At'*W*Ax];
   Ntt=[At'*W*At];
   
   
   A=[Ax At];
  
   
%     Nxx1=A'*W*A;
   N_xx=Nxx-Nxt*inv(Ntt)*Ntx;
   %Absolute vector
   nx=Ax'*W*dLf;
   nt=At'*W*dLf;
%    Reduced absolute vector
   n_x = nx -Nxt*inv(Ntt)*nt;
%   Create the G Constraint matrix
for i=1;j=1:2:24;
    G(i,j)= 1;
    i=2; j=2:2:24;
    G(i,j)=1;
    i=3;j=1:2:24;
    G(i,j)=Y;
    i=3;j=2:2:24;
    G(i,j)=-X;
    
end


%    Application of the G constrain Matrix
Gf=[G zeros(3,3)];
GM=[N_xx G';Gf];
%Solution for ?x
n_xf=[n_x;zeros(3,1)];
%From inv(GM) we will be able to get our Qxx 
delX=inv(GM)*n_xf;
% delT=inv(Ntt)*(nt-Ntx*delX);
%while delX>0.02
for j=1;
    i=1:2:24
    dX=delX(i);
    i=2:2:24;
    dY=delX(i);
end
%   disp('Value still too large')
%end


  X = X+dX;
  Y = Y+dY;
   
   dXf=[delX;zeros(2,1)];
   V=dLf-A*dXf;
  %end of the iterarative for loop
% end
   
   Nxx1=inv(GM);
   Qxx=Nxx1(1:24 ,1:24);
   
   
   %  aposteriori variance
     v = V'*W*V/(55-24);
     
% covariance matrix
     Exx = v * Qxx;
%      Exx_ = fliplr(Exx);
     E1 = diag(Exx);
     E2 = diag(Exx,1);
    
%  Variance of X and Y    
     for j=1;
         i=1:2:24;
          Varx =(E1(i));
          i1=1:2:23;
          Covarx = (E2(i1));
          i1=1:2:23;
         i=2:2:24;
          Covary = (E2(i1));
          Vary =(E1(i));
     end
% standard deviation of X and Y
 SDx = sqrt(Varx);
 SDy = sqrt(Vary);
 
a=sqrt((0.5*(Varx+Vary))+sqrt(0.25*(Varx-Vary)).^2+Covarx.^2);
b=sqrt((0.5*(Varx+Vary))-sqrt(0.25*(Varx-Vary)).^2+Covarx.^2);

        
%direction of the semimajor axis
     for i = 1:12;
         j=1;
        tanthita(i) = (2*Covarx(i))/(Varx(i)-Vary(i));
     end
        thita = (atan(tanthita)/2)';
        
     %plotting error ellipse
     for i = 1:12;
            x=X(i);
            y=Y(i);
            grid on;
            hold on;
            plot(x,y);
     end
    xlabel('X coordinates');
    ylabel('Y coordinates');
    title('ERROR ELLIPSE GROUP 1');
    
    for i = 1:49;
         p = X(Spd(i));
         g = X(Tpd(i));
         q = Y(Spd(i));
         r = Y(Tpd(i));
         x0 = [p g];
         y0 = [q r];
         plot(y0,x0); 
        
    end
    
    for i = 1:12;
        p = 1:1:12;
        t = -2*pi:pi/100:2*pi;
        x = X(i)+((a(i)*sin(t))*50);
        y = Y(i)+((b(i)*cos(t))*50);
        plot(y,x);
        text(Y(i),X(i),num2str(i));
    end

    %Global model test, we compute the test statistic
    
% Outlier detection by method of outliers
   
%    sigmaV = inv(W)-A*(inv(A'*W*A))*A';
   
%    sigmaVv = inv(W)-Ax*inv(Ax'*W*Ax)*Ax';
 
CSum= sum(X);
% SumXsq=sum((X(i)^2));

   
   
  
   
   