clear all;
%%%% %%%%%%%%%%%%% Author %%%%%%%%%% March 2019 %%%%%%%%%%%%%%%%%%
            %%% -ING. MWONGELA D. MATHINA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%create a file reader to read observations from excel file
format long;
Obsv=xlsread('ChemositDataV2.xlsx');
 Yobsv=Obsv(:,1);
 Xobsv=Obsv(:,2);
 Zobsv=Obsv(:,3);
 %%
 %Load Control Points from excel file
 Yctrl=Obsv(:,4);
 Xctrl=Obsv(:,5);
 Zctrl=Obsv(:,6);
 %%
 %%
 %   for loops that read values of X,Y and Z that are greater than zero    
for yy = 1:length(Yctrl); 
    
    if (Yctrl(yy)>0); 
        P(yy)=Yctrl(yy);
    end 
end
Yc=P';
for xx = 1:length(Xctrl); 
    
    if (Xctrl(xx)>0); 
        Q(xx)=Xctrl(xx);
    end 
end
Xc=Q';
for zz = 1:length(Zctrl); 
    
    if (Zctrl(zz)>0); 
        U(zz)=Zctrl(zz);
    end 
end
Zc=U';
%%
 %customized matrices to pick control and observations data
 rr=(1:1:33);
 R=rr'; %generate column matrix to pick observations
 bb=repmat(1:3:31,[3 1]);%generate a matrix to place A matrix elements in columns
 B=bb(:);
 G=1:3;
 gg=repmat(G,1,11); %generate column matrix to pick the control points 
 M=gg(:);   
 V=(1:3:97);%Place Ly elements in respective rows
 vv=V';
 S=(1:3:97);%Place Lx elements in respective rows
 ss=S';
 T=(1:3:97);%Place Lz elements in respective rows
 tt=T';
 K=(1:1:33);%Pick elements of the respective L-matrices(Ly,Lx and Lz) and place them in correct rows 
 kk=K';
 D=(1:3:97);%Place elements of the A-matrix in the resective rows
 dd=D';
 E=(1:1:33);%Placing adjusted coordinates in respective rows
 ee=E';
 F=repmat(1:1:11,[3 1]);%Pick adjusted coordinates
 ff=F(:);
 
 H=repmat(1,[11 1]);
 C=(1:1:11);
 cc=C(:);
 
 
 %%
 %Obtain baselines
 for i=1:33;
     m=M(i); %pick the control points
     r=R(i); %pick observations of the other points 
    dy=Yc(m)-Yobsv(r); %Obtain baselines aligned to Y direction
    dx=Xc(m)-Xobsv(r); %Obtain baselines aligned to X direction
    dz=Zc(m)-Zobsv(r); %Obtain baselines aligned to Z direction
    
    %baseline derived from control point? sum up with the control point 
    dyf(i,1)=-dy+Yc(m);
    dxf(i,1)=-dx+Xc(m);
    dzf(i,1)=-dz+Zc(m); 
    
    %Weight Matrix (1cm+2ppm.S)
    Wyc(i,1)=inv((0.01+((2/1000000)*abs(dy)))^2);
    Wxc(i,1)=inv((0.01+((2/1000000)*abs(dx)))^2);
    Wzc(i,1)=inv((0.01+((2/1000000)*abs(dz)))^2);
 end
%  
%%
%Extract the respective L-matrices
 Ly=dyf; %Y-direction L matrix
 Lx=dxf; %X-direction L matrix
 Lz=dzf; %Z-direction L matrix
 %The weight matrix
 Wy=Wyc;
 Wx=Wxc;
 Wz=Wzc;
 %%
 %for loop to organize and combine the distinct L-matrices 
for i=1:99;
   for i=1:33;
    k=kk(i);  %pick Ly elements
    v=vv(i);  %place Ly elements to resective rows
    i=v;
    L(i,1)=Ly(k);
    
    Wf(i,1)=Wyc(k);
   end
   for i=1:33;
    k=kk(i); %pick Lx elements
    s=ss(i); %place Lx elements to resective rows
    i=s;
    L(i+1,1)=Lx(k);
    
    Wf(i+1,1)=Wxc(k);
   end
   for i=1:33; 
    k=kk(i); %pick Lz elements
    t=tt(i); %place Lz elements to resective rows
    i=t;
    L(i+2,1)=Lz(k);
    
    Wf(i+2,1)=Wzc(k);
   end   
end
 %%
  %Weight Matrix design
  W=diag(Wf);
 %%
 %Design of the mighty A matrix
 for i=1:99;
        
        for i=1:33;
         b=B(i);    %place A matrix elements in the respective columns
         d=dd(i);   %place A matrix elements in the respective rows
         i=d;
         j=b;
         A(i,j)=1;
        end
        for i=1:33;
         b=B(i);    %place A matrix elements in the respective columns
         d=dd(i);   %place A matrix elements in the respective rows
         i=d;
         j=b;
         A(i+1,j+1)=1;
        end
        for i=1:33;
         b=B(i);    %place A matrix elements in the respective columns
         d=dd(i);   %place A matrix elements in the respective rows
         i=d;
         j=b;
         A(i+2,j+2)=1;
        end
       
 end
 %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
NN=A'*W*A; %%%%%Normal equation matrix %%%%%%%
Qxx=inv(NN);%%%%%%%% coffactor matrix %%%%%%
d1=A'*W*L;   %%%%% absolute vector %%%%%%%%
Soln=Qxx*d1; %%%%% coordinates %%%%%%%%
xlswrite('Soln', Soln)
%seperating individual coordinates
for j=1;
      i=1:3:33;
      dy1=Soln(i);
      i=2:3:33;
      dx1=Soln(i);
      i=3:3:33;
      dz1=Soln(i);
end
%% 
Y=dy1;
X=dx1;
Z=dz1;
%computing the residual vectors v;(v=L-Ax)
A_x=A*Soln;
v=L-A_x;

%%
% aposteriori variance sigma
sigma=(v'*v)/(99-33);
%%  
%Computing the standard deviation of the residuals
vmean= (sum(v))/99;
vsq= power(v-vmean,2);
Sd = sqrt((sum(vsq))/99);
%%
% perform Global Model Test GMM
 GMM=(v'*W*v);
%  Compute the covriance matrix
 Exx=sigma*Qxx;
%  Extracting covariance matrices associated with the 12 points
  E1 = diag(Exx);
  E2 = diag(Exx,1);
  
%  a for loop extracting the variances from the above matrix
 for j=1;
     i=1:3:33;
     sigx2=(E1(i));
     i1=1:3:32;
     Covarx = (E2(i1));
      i1=2:3:32;
     i=2:3:33;
     Covary = (E2(i1));
     sigy2= (E1(i));
 end
%%  
% obtaining elements of the error ellipses
 sigx=sqrt(sigx2);
 sigy=sqrt(sigy2);
 a=sqrt((0.5*(sigx2+sigy2))+sqrt((0.25*(sigx2-sigy2)).^2+Covarx.^2));
 b=sqrt((0.5*(sigx2+sigy2))-sqrt((0.25*(sigx2-sigy2)).^2+Covarx.^2));
 
 %  computation of the direction of semi-major axes
     for i = 1:11;
         j=1;
        tanthita(i) = (2*Covarx(i))/(sigx2(i)-sigy2(i));
     end
        thita = (atan(tanthita)/2)';
        
 %Plotting error ellipses
 for i = 1:11;
            x=X(i);
            y=Y(i);
            
            for i=1:3
              x1=Xctrl(i);  
              y1=Yctrl(i);  
            end
            
            grid on;
            hold on;
            plot(y,x);
            plot(y1,x1);
 end
    xlabel('X coordinates');
    ylabel('Y coordinates');
    title('ERROR ELLIPSES');
    
    %displaying the rays
 for i = 1:11;
         p = Xctrl(H(i));
         g = X(cc(i));
         q = Yctrl(H(i));
         r = Y(cc(i));
         x0 = [p g];
         y0 = [q r];
         plot(y0,x0);       
 end
    
 %     controlling the size of the ellipses
    for i = 1:11;
        p = 1:1:11;
        t = -2*pi:pi/50:2*pi;
        x = X(i)+((a(i)*sin(t))*100000);
        y = Y(i)+((b(i)*cos(t))*100000);
        plot(y,x);
        text(Y(i),X(i),num2str(i));
    end

 