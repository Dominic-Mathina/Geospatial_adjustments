clear all;
%%%% %%%%%%%%%%%%% Author %%%%%%%%%% March 2019 %%%%%%%%%%%%%%%%%%
            %%% -ING. MWONGELA D. MATHINA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%create a file reader to read observations from excel file
format long;
Obsv=xlsread('ChemositData.xlsx');
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
        S(zz)=Zctrl(zz);
    end 
end
Zc=S';
%%
 %customized matrices to pick control and observations data
  N=[1;2;3;1;2;3;1;2;3;4;5;6;4;5;6;4;5;6;7;8;9;7;8;9;7;8;9;10;11;12;10;11;12;10;11;12;13;14;15;13;14;15;13;14;15;16;17;18;16;17;18;16;17;18;19;20;21;19;20;21;19;20;21;22;23;24;22;23;24;22;23;24;25;26;27;25;26;27;25;26;27;28;29;30;28;29;30;28;29;30;31;32;33;31;32;33;31;32;33];
 bb=repmat(1:3:31,[9 1]);%generate a matrix to place A matrix elements in columns
 B=bb(:);
 G=1:3;
 gg=repmat(G,3,11);
 M=gg(:);
 hh=repmat((1:1:33),[3 1]);
 H=hh(:);
 ay=H';
  
 %Sequence numbers to place elements in rows
 V=(1:3:295);%Place Ly elements in respective rows
 vv=V';
 
 S=(1:3:295);%Place Lx elements in respective rows
 ss=S';
 T=(1:3:295);%Place Lz elements in respective rows
 tt=T';
 K=(1:1:99);%Pick elements of the respective L-matrices(Ly,Lx and Lz) and place them in correct rows 
 kk=K';
 D=(1:3:295);%Place elements of the A-matrix in the resective rows
 dd=D';
 E=(1:1:33);%Placing adjusted coordinates in respective rows
 ee=E';
 F=repmat(1:1:11,[3 1]);%Pick adjusted coordinates
 ff=F(:);
 
 
 %%
for i=1:3;
 %Obtain baselines
 for i=1:99;
     m=M(i); %pick the control points
     n=N(i); %pick observations of the other points 
    dy=Yc(m)-Yobsv(n); %Obtain baselines aligned to Y direction
    dx=Xc(m)-Xobsv(n); %Obtain baselines aligned to X direction
    dz=Zc(m)-Zobsv(n); %Obtain baselines aligned to Z direction
    
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
 
 Wy=Wyc;
 Wx=Wxc;
 Wz=Wzc;
 %%
 %for loop to organize and combine the distinct L-matrices 
for i=1:297;
   for i=1:99;
    k=kk(i);  %pick Ly elements
    v=vv(i);  %place Ly elements to resective rows
    i=v;
    L(i,1)=Ly(k);
    
    Wf(i,1)=Wyc(k);
   end
   for i=1:99;
    k=kk(i); %pick Lx elements
    s=ss(i); %place Lx elements to resective rows
    i=s;
    L(i+1,1)=Lx(k);
    
    Wf(i+1,1)=Wxc(k);
   end
   for i=1:99; 
    k=kk(i); %pick Lz elements
    t=tt(i); %place Lz elements to resective rows
    i=t;
    L(i+2,1)=Lz(k);
    
    Wf(i+2,1)=Wzc(k);
   end   
end
 %%
  %Weight Matrix
  W=diag(Wf);
 %%
 %Design of the mighty A matrix
 for i=1:297;
        
        for i=1:99;
         b=B(i);    %place A matrix elements in the respective columns
         d=dd(i);   %place A matrix elements in the respective rows
         i=d;
         j=b;
         A(i,j)=1;
        end
        for i=1:99;
         b=B(i);    %place A matrix elements in the respective columns
         d=dd(i);   %place A matrix elements in the respective rows
         i=d;
         j=b;
         A(i+1,j+1)=1;
        end
        for i=1:99;
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
for i=1:99;
    for i=1:33;
    i=ee(i);
    f=ff(i);
    Y(i,1)=dy1(f);
    end
    for i=1:33;
    i=ee(i);
    f=ff(i);
    X(i,1)=dx1(f);
    end
    for i=1:33'
    i=ee(i);
    f=ff(i);
    Z(i,1)=dz1(f);
    end
end
% Yobsv=Y;
% Xobsv=X;
% Zobsv=Z;
end 
%computing the residual vectors v;(v=L-Ax)
A_x=A*Soln;
v=L-A_x;

%%
% aposteriori variance sigma
sigma=(v'*W*v)/(297-33);
%%  
%Computing the standard deviation of the residuals
vmean= (sum(v))/297;
vsq= power(v-vmean,2);
Sd = sqrt((sum(vsq))/297);
%%
% perform Global Model Test GMM
 GMM=(v'*v);
%  Compute the covriance matrix
 Exx=sigma*Qxx;
%  Extracting covariance matrices associated with the 12 points
  E1 = diag(Exx);
  E2 = diag(Exx,1);

 
 
