
% image coordinates in principal point system [x-y mm]
[   4    12.472    -3.773;
   8    12.344     2.280;
  16     9.080    -7.899;
  41     7.983     7.573;
  55     4.772    -9.715;
 107    -0.834    10.514;
  123    -0.601    -9.906;
 175    -8.914     8.477;
 189    -5.751    -8.343;];
  
focallength=28.556;
% control points coordinates
           [ 4       5530.8       6576.1       3156.7
            8       6317.8       6539.4       2926.2
           16       4743.9       6553.2         2946
           41       6814.7       6466.8       2416.5
           55       4223.1         6485       2453.9
          107       7022.4       6361.7         1629
          123       3982.7       6374.3       1679.2
          175       6790.5       6262.5       855.51
          189       4187.4       6260.5       897.02];
% initial exterior orienation parameters
 w=90;phi=-33;k=92;xo=7700;yo=1817;zo=8900;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
photo=[   12.472    -3.773
    12.344     2.280;
    9.080    -7.899;
    7.983     7.573;
    4.772    -9.715;
   -0.834    10.514;
   -0.601    -9.906;
   -8.914     8.477;
   -5.751    -8.343];
xp=photo(:,1);yp=photo(:,2);delta=[1 1 1 1 1 1 1];

f=28.556;

ng=size(photo,1); 
 XYZ=[ 5530.8       3156.7       6576.1
       6317.8       2926.2       6539.4
       4743.9         2946       6553.2
       6814.7       2416.5       6466.8
       4223.1       2453.9         6485
       7022.4         1629       6361.7
       3982.7       1679.2       6374.3
       6790.5       855.51       6262.5
       4187.4       897.02       6260.5];
x=XYZ(:,1);y=XYZ(:,3);z=XYZ(:,2);
 
  
omega=(pi/180)*(90);phi=(pi/180)*(-33);kappa=(pi/180)*(92);xo=7700;yo=1817;zo=8900;
disp('INITIAL EXTERIOR ORIENTATION')
  wpk=[omega,phi,kappa,xo,yo,zo];
  [omega,phi,kappa,xo,yo,zo] 
  
%   [ Tx, Ty, Tz, w2, p2, k2 ]= Imageresection (XYZ,xp,yp,wpk,f );
  
  
 
 
 omega=wpk(1,1);phi=wpk(1,2);kappa=  wpk(1,3)  ;xo=wpk(1,4);yo=wpk(1,5);zo=wpk(1,6);
 xp=xx;yp=yy;delta=[1 1 1 1 1 1 1]; 
 
 
ng=size(XYZ,1); 
x=XYZ(:,1);y=XYZ(:,2); z=XYZ(:,3);
 
% _____________________________________________________
 ii=0;
 while max(abs(delta)) >.00001
     ii=ii+1;
%%%%%%%%%%%%%%%%%%%%% rotation matrix %%%%%%%%%%%%%%%%%%%%%
mw=[1 0 0;0 cos(omega) sin(omega) ;0 -sin(omega) cos(omega)];
mp=[cos(phi) 0 -sin(phi);0 1 0;sin(phi) 0 cos(phi)];
mk=[cos(kappa) sin(kappa) 0;-sin(kappa) cos(kappa) 0;0 0 1];
m=mk*mp*mw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% partial derivatives    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gg=ng*2;
for k=1:ng
    dx(k,1)=x(k,1)-xo;
    dy(k,1)=yo-y(k,1);
    dz(k,1)=z(k,1)-zo;
    q(k,1)=m(3,1)*(x(k,1)-xo)+m(3,2)*(z(k,1)-zo)+m(3,3)*(yo-y(k,1));
    r(k,1)=m(1,1)*(x(k,1)-xo)+m(1,2)*(z(k,1)-zo)+m(1,3)*(yo-y(k,1));
    s(k,1)=m(2,1)*(x(k,1)-xo)+m(2,2)*(z(k,1)-zo)+m(2,3)*(yo-y(k,1));
end
j=0;
for k=1:2:gg
    j=j+1;
    ff(k,1)=-(q(j,1)*xp(j,1)+r(j,1)*f)/q(j,1);
    ff(k+1,1)=-(q(j,1)*yp(j,1)+s(j,1)*f)/q(j,1);
end

j=0;
for k=1:2:gg
    j=j+1;
    b(k,1)=(xp(j,1)/q(j,1))*(-m(3,3)*dz(j,1)+m(3,2)*dy(j,1))+(f/q(j,1))*(-m(1,3)*dz(j,1)+m(1,2)*dy(j,1));
    b(k,2)=(xp(j,1)/q(j,1))*(dx(j,1)*cos(phi)+dz(j,1)*(sin(omega)*sin(phi))+dy(j,1)*(-sin(phi)*cos(omega)))+...
        (f/q(j,1))*(dx(j,1)*(-sin(phi)*cos(kappa))+dz(j,1)*(sin(omega)*cos(phi)*cos(kappa))+dy(j,1)*(-cos(omega)*cos(phi)*cos(kappa)));
    b(k,3)=(f/q(j,1))*(m(2,1)*dx(j,1)+m(2,2)*dz(j,1)+m(2,3)*dy(j,1));
    b(k,4)=-((xp(j,1)/q(j,1))*m(3,1)+(f/q(j,1))*m(1,1));
    b(k,5)=-((xp(j,1)/q(j,1))*m(3,2)+(f/q(j,1))*m(1,2));
    b(k,6)= ((xp(j,1)/q(j,1))*m(3,3)+(f/q(j,1))*m(1,3));
    b(k+1,1)=(yp(j,1)/q(j,1))*(-m(3,3)*dz(j,1)+m(3,2)*dy(j,1))+(f/q(j,1))*(-m(2,3)*dz(j,1)+m(2,2)*dy(j,1));
    b(k+1,2)=(yp(j,1)/q(j,1))*(dx(j,1)*cos(phi)+dz(j,1)*(sin(omega)*sin(phi))+dy(j,1)*(-sin(phi)*cos(omega)))+...
        (f/q(j,1))*(dx(j,1)*(sin(phi)*sin(kappa))+dz(j,1)*(-sin(omega)*cos(phi)*sin(kappa))+dy(j,1)*(cos(omega)*cos(phi)*sin(kappa)));
    b(k+1,3)=(f/q(j,1))*(-m(1,1)*dx(j,1)-m(1,2)*dz(j,1)-m(1,3)*dy(j,1));
    b(k+1,4)=-((yp(j,1)/q(j,1))*m(3,1)+(f/q(j,1))*m(2,1));
    b(k+1,5)=-((yp(j,1)/q(j,1))*m(3,2)+(f/q(j,1))*m(2,2));
    b(k+1,6)= ((yp(j,1)/q(j,1))*m(3,3)+(f/q(j,1))*m(2,3));

end
 format short g 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% Least Square %%%%%%%%%%%%%%%%%%%%5
  btb=inv(b'*b); 
  btf=b'*ff;
delta=btb*btf;v=b*delta-ff;D(:,ii)=delta;
omega=omega+delta(1,1);
phi  =phi  +delta(2,1);
kappa=kappa+delta(3,1);
xo=xo+delta(4,1);
yo=yo+delta(6,1);
zo=zo+delta(5,1);
end
sigm=sqrt(v'*v)/(size(b,1)-size(b,2));  
figure
subplot(1,2,1)
plot((1:size(D,2)), D(1,:),'-ro','Markerfacecolor','r','LineWidth',2),hold on
plot((1:size(D,2)), D(2,:),'-b^','Markerfacecolor','b','LineWidth',2),hold on
plot((1:size(D,2)), D(3,:),'-ms','Markerfacecolor','k','LineWidth',2),hold on
title(' LEAST SQUARE ADJUSTMENT - ANGLES CONVERGENCE') 
xlabel('NUMBER OF ITERATIONS')
ylabel('CORRECTIONS OF ORIENTATION ANGLES - RAD.')
legend('omega','phi','kappa')
grid on
axis tight
subplot(1,2,2)
plot((1:size(D,2)), D(4,:),'--M+','Markerfacecolor','m','LineWidth',2),hold on
plot((1:size(D,2)), D(5,:),'-ko','Markerfacecolor','c','LineWidth',2),hold on
plot((1:size(D,2)), D(6,:),'-b^','Markerfacecolor','r','LineWidth',2),hold on
title(' LEAST SQUARE ADJUSTMENT - CAMERA COORDINATES CONVERGENCE') 
xlabel('NUMBER OF ITERATIONS')
ylabel('CORRECTIONS OF CAMERA COORDINATES')
legend('Xo','Yo','Zo')
grid on
axis tight
w2=omega;
p2=phi;
k2=kappa;
Tx=xo;Ty=yo;Tz=zo; 
 
% disp('corrections for last iteration is:-');
% delta
disp('                                             ');
disp('EXTERIOR ORIENTATION PARAMETERS are:- ')
disp('*******************************************************');
disp(['adjusted omega (DEG.)=',num2str((180/pi)*omega)]);
disp(['adjusted phi   (DEG.)=',num2str((180/pi)*phi)]);
disp(['adjusted kappa (DEG.)=',num2str((180/pi)*kappa)]);
disp(['adjusted xo(m)=',num2str(xo)]);
disp(['adjusted yo(m)=',num2str(yo)]);
disp(['adjusted zo(m)=',num2str(zo)]);
disp('*******************************************************');
 
  

 