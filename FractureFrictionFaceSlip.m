clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 0: Bits you do not need to touch. Just leave these on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   


	%Script to run dislocations on a fault surface

	% This follows the same unit convention as Poly3d so ''assumes dimensionally consistent units for physical quantities with
	% dimensions of length or stress. Thus if you use kilometers for any quantity with
	% dimensions of length (e.g. a coordinate position), you must use kilometers for all
	% quantities with dimensions of length (e.g. displacement discontinuity).''
	% Andrew Lyell Thomas, Masters thesis, Stanford University. 

cmap = colormap_cpt('Ccool-warm'); %Loads a colourmap to be used for figures. This is Kenneth Morelands diverging colourmap.
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
%STEP 1: Import the fault surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  

 %%
r=3; %Circle rad
HalfLength=1; %Crack half length
Density=15;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Option C = User defined surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InnerDist=(1/(Density));
OuterEdgeSmpl=((2*pi)/InnerDist);

% OuterEdgeSmpl=round(rad2deg((2*pi)/(Density/2)));
OuterEdgeSmpl=round(OuterEdgeSmpl/2)*2; %Making sure its divible by two. 
Radius=1; 
X = linspace(-Radius,Radius,Density);
Y = linspace(-Radius,Radius,Density); 
[X,Y] = meshgrid(X,Y); 
%Inner circle of points
[Th,R]=cart2pol(X,Y); X(R>0.98)=[];Y(R>0.98)=[]; %Circle

[ xe,ye ] = CreateCircleXY( OuterEdgeSmpl,1,0 ); 
X=[X';xe]; Y=[Y';ye];
%Outer circle of points
[ xe,ye ] = CreateCircleXY( OuterEdgeSmpl/2,1+(2*InnerDist),0 ); 
%Rotating so we get equilateral tris. 
Theta=deg2rad(360/OuterEdgeSmpl)/4;
[Xrot,Yrot] = RotateObject2d(X,Y,Theta);
X=[X;xe]; Y=[Y;ye];
%Scaling so radius is one
X=X/(1+(2*InnerDist));
Y=Y/(1+(2*InnerDist));

[ Triangles,Points ] = MeshSurfaceXYPnts( X,Y );

%Rotating so it lies along X axis
[Points(:,2),Points(:,4)] = RotateObject2d(Points(:,2),Points(:,4),deg2rad(90));
[Points(:,2),Points(:,3)] = RotateObject2d(Points(:,2),Points(:,3),deg2rad(90));
% %Rotating around 90 so midpoint of fracture tip matches shear direction
% %(Makes less assymetric graphs)
% [Points(:,2),Points(:,4)] = RotateObject2d(Points(:,2),Points(:,4),deg2rad(90));


%% Cleaning mesh!
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
TR = triangulation(Triangles,Points(:,2:4)); 
[P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector]=CleanEdgeTris(MidPoint,P1,P2,P3,TR,FaceNormalVector);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2: Define Full/Halfspace and Elastic constants 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %Defining halfspace. 1 is on, 0 is off. The code will run much faster with this off as the calculations are simpler in a full space. 
	%Use this if the algorithm is throwing errors for vertex's above 0m.
    
halfspace = 0; 

	%Defining the height of the freesurface relative to the value of 0 of
	%your imported surfaces. Used when deformation is not at
	%sealevel. See the output MaxZ height if your half space clashes with
	%your points and adjust until this outputs as a negative value. 
freesurface_height = 0;

    %Defining elastic constants, All three are required to run this whole script.  See equations below if you need to calculate these 
	%from other constants. 
	%Note nu is needed in the displacement calculation.  
	
mu=120;%0.4;%4000 	%Shear Mod, mu or G. Relates shear stress to shear strain. 
nu = 0.25;     		%Poisson's ratio, Nu or V. Rubber 0.5, Cork 0, Rock 0.1-0.3;
[ K,E,lambda,nu,mu ] = ElasticConstantsCheck( nu,mu );

%Function to extract the three vertex's for each triangle and 
	%triangles defined by a combination of vertex Id's. This is used later by the algorithms 
    %The max Z value is output as text in the command window giving an idea
    %of minimum free surface height. 
Points=[Points(:,1),Points(:,2),Points(:,3),(Points(:,4)-freesurface_height)];
[MidPoint,FaceNormalVector] = MidPointCreate(Points,Triangles);
[P1,P2,P3] = CreateP1P2P3( Triangles,Points );
	%Fictitious disp flag, if set to one then this forces the triangles with this to 0 displacement. 
Fdisp= zeros(size(Triangles(:,1))); 

close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 3: Define how you want to calculate slip - Constant, Remote, Remote
%no opening, tractions on the elements or user imported slip value for each face. 

%Stress tensor outputs from this part of the toolbox i.e. Sxx Syy etc are simply the user defined remote stress
%components (boundary conditions). These are used later in the stress calculations. If strain is input these are converted to stress 
%using Hooke's law. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[Sxx,Syy,Szz,Sxy,Sxz,Syz,Tnn,Tss,Tds,Mu,Sf,strain ] = CreateBlankVars;
	%%%%%%%%%%%%%%
    %StrainInput        If halfspace = 1, Z stress components are removed
	%%%%%%%%%%%%%%
    
strain=0;                   %Put to 1 to define the stresses defined in 'stress input' as strain values
    

	%%%%%%%%%%%%%% 
    %StressInput
	%%%%%%%%%%%%%%
	
ne=numel(MidPoint(:,3));
Sxx = -20.15; 					
Syy = -8.23;	
Szz = Sxx;
Sxy = -7.0;     		
Sxz = 0;
Syz = 0; 

Mu  =ones((ne),1)*0.6;    %Coefficient of friction
Sf  = ones(size(Mu))*0.0;  

Tnn=0;
Tds=0;
Tss=0;

Option='C'; %slip from uniform remote stress with friction
%Sxx=0; Syy=0; Szz=0; Sxy=0; Sxz=0; Syz=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that does the work:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Dss,Dds,Dn,Sxx,Syy,Szz,...
Sxy,Sxz,Syz]=SlipCalculator3d(MidPoint,Sxx,Syy,Szz,Sxy,Sxz,Syz,...
Tnn,Tss,Tds,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option);

%Reset some vars:
FlatPlaneNorm=[0,-1,0];
[ Tnn,Tds,Tss ] = CalculateNormalAndShearTractions3d( FlatPlaneNorm,Sxx,Syy,Szz,Sxy,Sxz,Syz );
TsDr=Tss-abs(Tnn*Mu(1));
%Analytical
VolPred=(8*(nu-1)*TsDr(1)*Radius^3)/(3*(nu - 2)*mu);
AreaOfCirc=pi*Radius^2;
Vol1AnResult=VolPred/AreaOfCirc;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Drawing figures of slip distribution on the crack.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draws 3 figures of the slip distribution on the surfaces
PlotSlipDistribution3d(Triangles,Points,cmap2,Dss,Dds,Dn)


%Numerical
[ Area,HPerim ] = AreaOfTriangle3d( P1(:,1),P1(:,2),P1(:,3),P2(:,1),P2(:,2),P2(:,3),P3(:,1),P3(:,2),P3(:,3) );
[ StrikeSlipDisp,Area ] = RowVecToCol( Dss,Area);
Vol=(sum(Area.*Dss)/2);
Vol1NumResult=Vol/sum(Area);
RevSlip=any(Dss<0);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP AA: Analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%An Eq
L=0:0.001:1;
a=1;
DnAn=((4*(1-nu)*a*Tnn)/(pi*mu))*sqrt(1-(L.^2/a^2));
DsAn=((8*(1-nu)*a*TsDr)/(pi*(2-nu)*mu))*sqrt(1-((L.^2)/(a^2)));

MxAn=max(DsAn);
%Normalising An for graphs
DsAn=DsAn./MxAn;
a=Radius;
%Normalising Num for graphs
Dss=(Dss)./MxAn;
[az,el,rdata]=cart2sph(MidPoint(:,1),MidPoint(:,2),MidPoint(:,3));
[Th,~]=cart2pol(MidPoint(:,1),MidPoint(:,3));

%Drawing the figure
pp=figure;
subplot(3,4,[1 2 3  5 6 7  9 10 11 ],'Position',[0.13,0.11,0.5276,0.815])
hold on
%scatter(R,TotalShearing,'.');
%Plotting analytical profiles
h1=plot(L,DsAn,'b','LineWidth',2.5);
%Plotting interpolated numerical disps every 20th halflength
scatter(rdata,Dss,24,'k','filled');
x1=max(rdata)-0.02;
x2=max(rdata)+0.02;
y1=min(Dss)-0.07;
y2=min(Dss)+0.07;
x = [x1, x2, x2, x1, x1];
y = [y1, y1, y2, y2, y1];
plot(x, y, 'k');
%Adding titles etc
title('Benchmarking slip distributions');
xlabel('\slL')
ylabel({'D_{\sls} / D_{\sls,p} _m_a_x'}) 
grid on
legend('show')
legend('Penny-shaped crack D_{\sls,p}','Numerical result')


subplot(3,4,[4 8])
plot(L,DsAn,'b','LineWidth',2.5);
hold on
grid on
%Plotting interpolated numerical disps every 20th halflength
scatter(rdata,Dss,24,'k','filled');
xlim([x1 x2])
ylim([y1 y2])
yticks([y1 (y2+y1)/2 y2])

%% Compute error
DnAn=((4*(1-nu)*a*Tnn)/(pi*mu))*sqrt(1-(rdata.^2/a^2));
DsAn=((8*(1-nu)*a*TsDr)/(pi*(2-nu)*mu))*sqrt(1-((rdata.^2)/(a^2)));

DsAn=DsAn./MxAn; %normalise

%Grab edge points
Edges=round(rdata,2)==round(max(rdata),2);
Error=((100./DsAn(Edges)).*Dss(Edges))-100;
subplot(3,4,[12])
scatter(rad2deg(abs(Th(Edges)-pi)),Error,24,'k','filled')
xlabel('\theta^{\circ}')
xlim([0 360])
xticks([0 360/2 360])
%ylabel('Error Perc')
grid on

%Making the plot nicer
titlesz=20;
fntsz=18;
ChangeFontSizes(fntsz,titlesz);
