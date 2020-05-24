
clear;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               space photo resection by collinearity equations         %%
%% Author:  Bashar Alsadik  - University of Twente  - ITC - 2010         %%
%% Based on curriculum book of "Analytical photogrammetry" - 2009        %%
%% at Baghdad University - College of engineering  - Surveying Dept.     %%
%% input:                                                                %%
%% - XYZ: Reference point coordinates                                    %%
%% - image coordinates [mm] :xx,yy [in p.p. system]                      %%
%% - f: focal length [mm]                                                %%
%% - wpk: initial exterior orienattion parameters                        %%
%% output:                                                               %%
%% - adjusted exterior orienattion parameters                            %%
%% Note: reference points should be at least (3 points)                  %%
%% PLEASE REFER TO THE CODE AUTHOR WHEN USING IT                         %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% image coordinates in principal point system [x-y mm]
%    4    12.472    -3.773;
%    8    12.344     2.280;
%   16     9.080    -7.899;
%   41     7.983     7.573;
%   55     4.772    -9.715;
%  107    -0.834    10.514;
%   123    -0.601    -9.906;
%  175    -8.914     8.477;
%  189    -5.751    -8.343;
  
% focal length=28.556 mm
% control points coordinates
%             4       5530.8       6576.1       3156.7
%             8       6317.8       6539.4       2926.2
%            16       4743.9       6553.2         2946
%            41       6814.7       6466.8       2416.5
%            55       4223.1         6485       2453.9
%           107       7022.4       6361.7         1629
%           123       3982.7       6374.3       1679.2
%           175       6790.5       6262.5       855.51
%           189       4187.4       6260.5       897.02
% initial exterior orienation parameters
%  w=90;phi=-33;k=92;xo=7700;yo=1817;zo=8900;
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
  
  [ Tx, Ty, Tz, w2, p2, k2 ]= Imageresection (XYZ,xp,yp,wpk,f );
  

 