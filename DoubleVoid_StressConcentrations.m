%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you dont need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   

clear;close all
	
	%Script to run dislocations on a fault surface

	% This follows the same unit convention as Poly3d so ''assumes dimensionally consistent units for physical quantities with
	% dimensions of length or stress. Thus if you use kilometers for any quantity with
	% dimensions of length (e.g. a coordinate position), you must use kilometers for all
	% quantities with dimensions of length (e.g. displacement discontinuity).''
	% Andrew Lyell Thomas, Masters thesis, Stanford University. 

cmap = colormap_cpt('Ccool-warm'); %Loads a colourmap to be used for figures. This is Kenneth Morelands divirging colourmap.
cmap2 = colormap_cpt('Ccool-warm2');

	%Conventions: 
		%See figure 6.14 in David Pollards book for tensors. 
		%Geological convention is that a confining pressure is a positive stress. Geological convention is used throughout this script. 
		%For a normal stress: positive if it produces compression in the material.
		%For a shear stress: negative if, when acting on the + face identified by the first index, it points
		%in the + direction identified by the second index. Example: Exy is - if on the +x face it points
		%in the +y direction.
		%For example a positive value in shear stress Sxy will cause RightLateral movement.
		
		%Strain uses the same convention as engineering stress convention where extension is a positve value. When entering strain values remember this. 
		%For a normal stress: positive (negative) if it produces tension (compression) in the material.
		%For a shear stress: positive if, when acting on the + face identified by the first index, it points
		%in the + direction identified by the second index. Example: Exy is + if on the +x face it points
		%in the +y direction.
		%For example a positive value in shear strain Exy will cause Leftlateral movement.
		
		
%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Loading Gocad ascii data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'DefaultFigureVisible','on');

Xmv=1.1; %works for all voids
numOpoints=numel(Xmv);
Ymv=zeros(size(Xmv));
Zmv=zeros(size(Xmv));

i=1;
p=1;
q=1;

xsclR=1;
ysclR=1;
zsclR=1;

 string='SphereFix4.ts';
 [ PointsF1,TrianglesF1 ]= GoCadAsciiReader( string ); 
  PointsF1(:,2:4)=PointsF1(:,2:4)*3;
  
string='SphereUniformDistribution_Eq2000FacesSides.ts';
[ Points1,Triangles1 ] = GoCadAsciiReader( string );

  %For later (defining location of obs pnts)
Sizer=numel(Triangles1(:,1));
  
%want rad as 0.5  
PointsF1(:,2)=PointsF1(:,2)*0.5;
PointsF1(:,3)=PointsF1(:,3)*0.5;
PointsF1(:,4)=PointsF1(:,4)*0.5;
Points1(:,2)=Points1(:,2)*0.5;
Points1(:,3)=Points1(:,3)*0.5;
Points1(:,4)=Points1(:,4)*0.5;

%scaling
PointsF1(:,2)=PointsF1(:,2)*xsclR;
PointsF1(:,3)=PointsF1(:,3)*ysclR;
PointsF1(:,4)=PointsF1(:,4)*zsclR;
Points1(:,2)=Points1(:,2)*xsclR;
Points1(:,3)=Points1(:,3)*ysclR;
Points1(:,4)=Points1(:,4)*zsclR;

figure; scatter3(Xmv(:),Ymv(:),Zmv(:),15)
LstMv=[Xmv(:),Ymv(:),Zmv(:)];    
  

%Movement
xmv=LstMv(i,1);
ymv=LstMv(i,2);
zmv=LstMv(i,3);

rx=xsclR/2; %the horizontal radii
ry=ysclR/2; %the vertical radii
rz=zsclR/2; %the vertical radii
%Tolerance, stress are not calculated closer than this to the edges and
%intersection flag is raised if they are closer than this.
tol_dist=0.1; 
zsclt=zsclR+tol_dist;
xsclt=xsclR+tol_dist;

%CHECK TO SEE IF THE TWO ELLIPSES WILL INTERSECT
%angle
theta=atan2(zmv,xmv);

angle_in_deg=rad2deg(theta);

%This length is the length of the line at the specified angle to the edge
%of the ellipse
RadiusAtAngle = sqrt(1/((cos(theta)/xsclt)^2+(sin(theta)/zsclt)^2));
RadiusAtAngle=RadiusAtAngle/2;

MidPointDistances = ((zmv)^2 + (xmv)^2)^0.5;

%Creating 2ndSphere
Points2=Points1;
Points2(:,2)=Points2(:,2)+xmv;%move
Points2(:,3)=Points2(:,3)+ymv;
Points2(:,4)=Points2(:,4)+zmv;

%Appending the spheres
BoundaryFlag=zeros(size(Triangles1(:,1)));
[Points,Triangles,BoundaryFlag] = DataAppender3d( Points1,Points2,Triangles1,Triangles1,BoundaryFlag,0 );

%Creating 2nd fixed pnts
PointsF2=PointsF1;
PointsF2(:,2)=PointsF2(:,2)+xmv;%move
PointsF2(:,3)=PointsF2(:,3)+ymv;
PointsF2(:,4)=PointsF2(:,4)+zmv;

%Appending the fixed points
[PointsF,TrianglesF,~] = DataAppender3d( PointsF1,PointsF2,TrianglesF1,TrianglesF1); clear dud %not needed


%Appending all
[Points,Triangles,Fdisp] = DataAppender3d( Points,PointsF,Triangles,TrianglesF,BoundaryFlag,1 );

 %the half length of each axis of the ellipsoid
 a=rx; 
 b=ry; 
 c=rz; 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calclations are simpler in a full space. 
    %Use this if the algorithm is throwing errors for vertex's above 0m.

halfspace = 0; 

    %Defining the height of the freesurface relative to the value of 0 of
    %your imported surfaces. Used when deformation is not at
    %sealevel. See the output MaxZ height if your half space clashes with
    %your points and adjust until this outputs as a negative value. 
freesurface_height =0;

    %Defining elastic constants, All three are required to run this whole script.  See equations below if you need to calculate these 
    %from other constants. 
    %Note nu is needed in the displacement calculation.  

mu=150000; %1%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

    %Function to extract the three vertex's for each triangle and 
    %triangles defined by a combination of vertex Id's. This is used later by the algorithms 
    %The max Z value is output as text in the command window giving an idea
    %of minimum free surface height. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
clear MidPoint FaceNormalVector
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
close
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 


%FOR DRAWING THE LOOP
hold on
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);axis equal;
%Forcing this to draw
drawnow
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hookes law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,Mu,Sf,strain,SecondSurface ] = CreateBlankVars;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option B = Run Influence code to see how the fault reacts to a remote
    %stress defined by the user. Choose to define in stress or strain. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%
    %StrainInput        If halfspace = 1, Z stress components are removed
    %%%%%%%%%%%%%%

strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values

    %%%%%%%%%%%%%% 
    %StressInput
    %%%%%%%%%%%%%%


Sxx = 0;                    %Positive is extension
Syy = 0; 
Szz = -1;
Sxy = 0;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz = 0;
Syz = 0; 
Option='B'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Going2SlipCalc')

if ~strcmp(Option,'A')    
    [Dss,Dds,Dn,Sxx,Syy,Szz,Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,...
     Syy,Szz,Sxy,Sxz,Syz,Tn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,...
     FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points);
end

disp('OutOfSlipCalc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing any fixed triangles
if any(Fdisp)==1
[Triangles,FaceNormalVector,MidPoint,P1,P2,P3]...
    = RemovingFixedEls3d(Triangles,FaceNormalVector,MidPoint,P1,P2,P3,Fdisp);
end

PlotSlipDistribution3d(Triangles,Points,cmap2,Dn,Dds,Dss)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %STEP 4: Define dispersed XYZ oberservation points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Option A = Simply define flat observation plane of points. XYZ with defined step size 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % %%%%Using midpoints%%%%
X=MidPoint(:,1);
Y=MidPoint(:,2);
Z=MidPoint(:,3);

%Drawing figure of surface and the observation points
figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
hold on
scatter3(X(:),Y(:),Z(:),5,'k')  %Showing surface and obs points
xlabel('x'); ylabel('y');zlabel('z'); axis('equal'); title('Surface and Obs Points');

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %STEP 5: Calculate Stresses  on dispersed XYZ
    %oberservation points. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Option A = Calculate Stresses and Strains
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TotalStress,Stress and StressChange are Xx12 arrays of strain tensors and stress tensors. XX YY ZZ XY XZ YZ
        %Strain is the first 6 coloumns and stress the last 6. Stress is the regional stress 
        %StressChange is the stress on the points from the event.
        %TotalStress is the driving stress and stress change from the event added. 

Sxx = 0;                    %Positive is extension
Syy = 0; 
Szz = -1;
Sxy = 0;        			%Positive = left lat when frac strikes NS and right lat when EW. See Pollard Fletcher Diagram - 6.13
Sxz = 0;% 0.05;
Syz = 0; 

[StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]=...
CalculateStressOnSurroundingPoints3d(Dss,Dds,Dn,mu,lambda,X,Y,Z,Sxx,...
Syy,Szz,Sxy,Sxz,Syz,P1,P2,P3,halfspace,nu);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %STEP 7: VISULISATION AND ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sxx,Syy,Szz,Sxy,Sxz,Syz ] = ExtractCols( StressTTotal );
[Exx,Eyy,Ezz,Exy,Exz,Eyz ] = ExtractCols( StrainTTotal );

%% LINES TO REMOVE

[ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( FaceNormalVector );
[ Srrnum,Sppnum,Sttnum,P12,P13,P23 ] = StressTensorTransformation3d(Sxx,Syy,Szz,Sxy,Sxz,Syz,FaceNormalVector,DipSlipCosine,StrikeSlipCosine);

figure;
trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceColor', [0.5 0 0.9 ]);axis equal;
hold on
quiver3(X,Y,Z,FaceNormalVector(:,1),FaceNormalVector(:,2),FaceNormalVector(:,3));
quiver3(X,Y,Z,StrikeSlipCosine(:,1),StrikeSlipCosine(:,2),StrikeSlipCosine(:,3));
quiver3(X,Y,Z,DipSlipCosine(:,1),DipSlipCosine(:,2),DipSlipCosine(:,3));

RadStr=Srrnum.*2;
PhiStr=Sppnum.*2;
ThStr=Sttnum.*2;
PlotSlipDistribution3d(Triangles,Points,cmap2,RadStr,PhiStr,ThStr)

[S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Srrnum*2,Sppnum*2,Sttnum*2,P12*2,P13*2,P23*2);
PlotSlipDistribution3d(Triangles,Points,cmap2,S1)

%Only works for one sampling this!
Voidtop=acosd(abs(FaceNormalVector(:,3)))<(1*zsclR);
Voidside=acosd(abs(FaceNormalVector(:,1)))<(1*xsclR) | acosd(abs(FaceNormalVector(:,2)))<(1*ysclR);

Sizze=numel(Triangles1)/3;
%sphere 1 pointing along X
PointingAlongX=round(acosd(FaceNormalVector(:,1)),9)==0;
PointingAlongX(Sizze+1:end)=0;
%sphere 2 pointing along -X
PointingAlongX2=round(acosd(FaceNormalVector(:,1)),9)==180;
PointingAlongX2(1:Sizze)=0;
Voidside=(PointingAlongX+PointingAlongX2)==1;



if sum(Voidtop)==0 || sum(Voidside)==0
    disp('using a bad sampling ')
    Smls=acosd(abs(FaceNormalVector(:,3)))-(1*zsclR);
    [~,Indx] = sort(Smls(:));
    VoidtopIdx = Indx(1:4);
    Voidside(VoidtopIdx)=1;
    Voidtop=Voidside;
    figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),double(Voidtop));
   
    Voidside=Voidside*0; %reset
    Smls=acosd(abs(FaceNormalVector(:,1)))-(1*xsclR); 
    Smls2=acosd(abs(FaceNormalVector(:,2)))-(1*ysclR);
    [~,Indx] = sort(Smls(:));
    VoidsdIdx = Indx(1:4);
    [~,Indx2] = sort(Smls2(:));
    VoidsdIdx = [VoidsdIdx;Indx2(1:4)];
    Voidside(VoidsdIdx)=1;
    Voidside=logical(Voidside);
    figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),double(Voidside));
else
    figure;trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),double(Voidside));
end


%most extensional
SttMax=max(max([Sttnum(Voidtop);Sppnum(Voidtop)]*2)); %not sure on if its spp or stt due to flat tri
%most compressional
SppMin=min(min(Sppnum(Voidside)*2)); 
%most compressional intermediate
IntMin=min(min(Sttnum(Voidside)*2));
%most extensional intermediate
IntMax=max(max(Sttnum(Voidside)*2));

disp('Max Tension')
disp(SttMax)
disp('Max Compression')
disp(SppMin)

set(0,'DefaultFigureVisible','on');



