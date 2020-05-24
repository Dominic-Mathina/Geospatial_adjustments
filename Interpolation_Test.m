clear all;
%%%% %%%%%%%%%%%%% Author %%%%%%%%%% March 2019 %%%%%%%%%%%%%%%%%%
            %%% -ING. MWONGELA D. MATHINA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%create a file reader to read observations from excel file
format long;
Obsv=xlsread('kbkData.xlsx');
 Bo=Obsv(:,1)*10;
 X=Obsv(:,2)/10000;
 Y=Obsv(:,3)/100000;
 Z=Obsv(:,4);
 
 Obsv1=xlsread('kbktst.xlsx');
 Xu=Obsv1(:,1);
 Yu=Obsv1(:,2); 
 %%
 %%
A = [Bo X Y ];

N=A'*A;

Qxx = inv(N);

d = A'*Z;

Soln = Qxx * d;
%Extract Values%
bo = Soln(1,1)*10;
b1 = Soln(2,1)/10000;
b2 = Soln(3,1)/100000;
%%

for i = 1:length(Xu)
    Zxy(i,1) = (bo) + (b1*Xu(i)) + (b2*Yu(i))
end

%%xlswrite('Zxy', Zxy)
