%%% pcolm - a simplified version of pcol.f code (Peric, 1997)
%%%
%%% This version is writen for matlab.
%%%
%%% Gabriel Usera - Feb/2007
%%%
%%% This version is setup for Lid-Driven Cavity flow.
%%%
%%% To set up other flows B.C. require modification.
%%%
%%% Typical user input is covered in lines 43 to 90,
%%% including all numerical parameters, fluid properties
%%% and grid specification (lines 89 and 90).
%%%
%%% Boundary Conditions are specified in lines :
%%% 126     - Initial B.C. setup for Velocity
%%% 288-332 - Complete B.C. for Velocity
%%% 431-435 - B.C. for mass balance/pressure equation
%%%
%%% Please keep this version 'as it is' and rename
%%% the versions that you modify ('poclm9538.m' or
%%% something like that).
%%%
%%% To run :
%%% >> [X,Y,XC,YC,FX,FY,U,V,P,F1,F2]=pcolm;
%%%
%%% To plot :
%%% >> plotm;
%%%
%%%
function [X,Y,XC,YC,FX,FY,U,V,P,F1,F2]=pcolm;

global LTIME LTEST                              % Logic variables, Steady/Transient and Testing
global IPR JPR IU IV IP                         % (IPR,JPR) Pressure reference. IU,IV,IP Equation tags.
global NI NJ NIM NJM NIJ LI                     % Domain dimensions and indexing vector LI=(I-1)*NJ
global X Y XC YC FX FY                          % Domain coordinates and interpolation factors
global DENSIT VISC GDS DT DTR SMALL ALFA ULID   % Fluid properties, Blending factor, Time step,...
global U V P PP F1 F2 DPX DPY U0 V0             % Variables: velocities, pressure, mass flux, pressure gradient...
global AE AW AS AN AP APR SU SV                 % Coefficient matrices
global SOR URF NSW RESOR                        % Iteration control parameters

%%% INPUT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ITIM=0;       % Set iteration counter to 0
TIME=0.;      % Set initial time to 0.
%
LTEST=0;      % Logic Flag : Testing ? LTEST=1, otherwise LTEST=0.
LTIME=0;      % Logic Flag : Stationary (0) or Transient (1) ?
%
MAXIT=500;    % Maximum number of outer iterations in each time step
IMON=5;       % Index (I) for monitoring location
JMON=5;       % Index (J) for monitoring location
IPR=2;        % Index (I) for pressure reference point
JPR=2;        % Index (J) for pressure reference point
SORMAX=1e-4;  % Residual level for stopping outer iterations
SLARGE=1e+3;  % Residual level for divergence
ALFA=0.92;    % Parameter for linear solver SIPSOL
%
DENSIT=1.0;   % Fluid density
VISC=  1.e-2; % Fluid viscosity
ULID=  1.;    % Lid velcity for lid driven cavity
%
UIN=0.;       % Initial value for U velocity (usually 0.)
VIN=0.;       % Initial value for V velocity (usually 0.)
PIN=0.;       % Initial valur for Pressure   (usually 0.)
%
ITST=1;       % Number of time steps (1 for steady flow, LTIME=0)
DT=1.00;      % Time step in seconds (meaningless for steady flow)
%
IU=1;         % Tag for U equation
IV=2;         % Tag for V equation
IP=3;         % Tag for P equation
%
URF(IU)=0.8;  % Under relaxation parameter for U
URF(IV)=0.8;  % Under relaxation parameter for V
URF(IP)=0.2;  % Under relaxation parameter for P
%
SOR(IU)=0.2;  % Ratio of residual reduction for linear solver : U
SOR(IV)=0.2;  % Ratio of residual reduction for linear solver : V
SOR(IP)=0.2;  % Ratio of residual reduction for linear solver : P
%
NSW(IU)=1;    % Maximum number of linear solver iterations : U
NSW(IV)=1;    % Maximum number of linear solver iterations : V
NSW(IP)=6;    % Maximum number of linear solver iterations : P
%
GDS=1.0;      % UDS - CDS Blending for U,V
%

%%% GRID DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[0:.025:1,1]'; % X - cell volume coordinates (only tested uniform spacing)
Y=[0:.025:1,1]'; % Y - cell volume coordinates (only tested uniform spacing)
%
NI=length(X);  % Number of points in X ( + 1 )
NJ=length(Y);  % Number of points in Y ( + 1 )
NIM=NI-1;      % NI - 1
NJM=NJ-1;      % NJ - 1
NIJ=NI*NJ;     % Total number of points 

%%% Indexing vector for k=(i-1)*NJ+j indexing style
LI=([1:NI]-1)*NJ;   % 

%%% X cell center coordinate
XC=X; XC(2:NIM)=0.5*(X(1:NIM-1)+X(2:NIM));
%%% Y cell center coordinate
YC=Y; YC(2:NJM)=0.5*(Y(1:NJM-1)+Y(2:NJM));

%%% Interpolation factors (in uniform grid FX=FY=0.5 except at boundaries)
FX(1:NIM)=(X(1:NIM)-XC(1:NIM))./(XC(2:NI)-XC(1:NIM)); FX(NI)=0.;
FY(1:NJM)=(Y(1:NJM)-YC(1:NJM))./(YC(2:NJ)-YC(1:NJM)); FY(NJ)=0.;

%%% SET SOME CONTROL VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SMALL=1.e-15; % Small number
GREAT=1.e+15; % Big number
DTR=1./DT;    % Reciprocal of time step

RESOR=zeros(3,1); % Initialize array to store residuals

%%% Monitoring location
IJMON=LI(IMON)+JMON;

%%% INITIALIZE FIELDS and ARRAYS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U,V,P,PP,U0,V0,F1,F2,DPX,DPY]=deal(zeros(NIJ,1));
[AP,AE,AW,AN,AS,APR,SU,SV]=deal(zeros(NIJ,1));

%%% FIXED BOUNDARY AND INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% North wall velocity
for I=2:NIM, U(LI(I)+NJ)=ULID; end;
%%% Intial conditions
for I=2:NIM,
    for IJ=LI(I)+2:LI(I)+NJM,
        U(IJ)=UIN;
        V(IJ)=VIN;
        P(IJ)=PIN;
        U0(IJ)=UIN;
        V0(IJ)=VIN;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ITIMS=ITIM+1;    % Update iteration counter to first iteration
ITIME=ITIM+ITST; % Final number of iterations
%
for ITIM=ITIMS:ITIME,           % Set time loop
    TIME=TIME+DT;               % Update time
    %
    if LTIME, U0=U; V0=V; end;  % Shift solutions in time
    %
    %%% HEADDING FOR THIS TIME STEP
    disp([' ']);
    disp(['**************************************************************']);
    disp(['TIME=',num2str(TIME,'%0.2E%')]);
    disp([' ']);
    disp(['IT.--RES(U)----RES(V)----RES(P)-----UMON-----VMON-----PMON----']);
    %
    %
    %%% OUTER ITERATIONS LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ITER=1:MAXIT
        CALCUV;         % Call CALCUV routine to update U and V
        CALCP;          % Call CALCP routine to update P,F1,F2 and U,V again
        %
        % Display information about residuals and monitor point values
        disp([num2str(ITER),'  ',num2str(RESOR(IU),'%0.2E'),' ',num2str(RESOR(IV),'%0.2e'),...
                             ' ',num2str(RESOR(IP),'%0.2e'),' ',num2str(U(IJMON) ,'%0.2e'),...
                             ' ',num2str(V(IJMON) ,'%0.2e'),' ',num2str(P(IJMON) ,'%0.2e')]);
        %
        % Check convergence
        SOURCE=max(RESOR([IU,IV,IP]));
        if SOURCE>SLARGE, break; end;
        if SOURCE<SORMAX, break; end;
    end;
    %
    % Return/Break if DIVERGING
    if SOURCE>SLARGE, disp('DIVERGING'); return; end;
end;
%
disp([' ']);
disp(['CALCULATION FINISHED - SEE RESULTS IN [X,Y,XC,YC,FX,FY,U,V,P,F1,F2]']);
%

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CALCUV

global LTIME LTEST
global IPR JPR IU IV IP
global NI NJ NIM NJM NIJ LI
global X Y XC YC FX FY
global DENSIT VISC GDS DT DTR SMALL ALFA ULID
global U V P PP F1 F2 DPX DPY U0 V0
global AE AW AS AN AP APR SU SV
global SOR URF NSW RESOR

%%% INITIALIZE COEFFICIENTS AND SOURCE TERMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AP,AE,AW,AN,AS,SU,SV,APR]=deal(zeros(NIJ,1));

%%% FLUXES THROUGH INTERNAL EAST CV FACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=2:NIM-1,
    FXE=FX(I);             % Interpolation Factor
    FXP=1.-FXE;            %
    DXPE=XC(I+1)-XC(I);    % Distance P->E
    %
    for J=2:NJM,
        IJ=LI(I)+J;        % Index IJ to P
        IJE=IJ+NJ;         % Index IJE to E
        %
        S=Y(J)-Y(J-1);     % Cell face 'area'
        D=VISC*S/DXPE;     % Coefficient for diffusive flux
        %
        CE=min(F1(IJ),0.); % Mass fluxes, upwind from E
        CP=max(F1(IJ),0.); % Mass fluxes, upwind from P
        %
        FUUDS=CP*U(IJ)+CE*U(IJE);            % Explicit convective fluxes UDS, U
        FVUDS=CP*V(IJ)+CE*V(IJE);            % Explicit convective fluxes UDS, V
        FUCDS=F1(IJ)*(U(IJE)*FXE+U(IJ)*FXP); % Explicit convective fluxes CDS, U
        FVCDS=F1(IJ)*(V(IJE)*FXE+V(IJ)*FXP); % Explicit convective fluxes CDS, V       
        %
        AE(IJ )=+CE-D;                       % Coefficient for E in P
        AW(IJE)=-CP-D;                       % Coefficient for W in E
        %
        SU(IJ )=SU(IJ )+GDS*(FUUDS-FUCDS);   % Source term for P, U
        SU(IJE)=SU(IJE)-GDS*(FUUDS-FUCDS);   % Source term for E, U
        SV(IJ )=SV(IJ )+GDS*(FVUDS-FVCDS);   % Source term for P, V
        SV(IJE)=SV(IJE)-GDS*(FVUDS-FVCDS);   % Source term for E, V
    end;
end;

%%% FLUXES THROUGH INTERNAL NORTH CV FACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for J=2:NJM-1,
    FYN=FY(J);             % Interpolation Factor
    FYP=1.-FYN;            %
    DYPN=YC(J+1)-YC(J);    % Distance P->N
    %
    for I=2:NIM,
        IJ=LI(I)+J;        % Index IJ to P
        IJN=IJ+1;          % Index IJN to N
        %
        S=X(I)-X(I-1);     % Cell face 'area'
        D=VISC*S/DYPN;     % Coefficient for diffusive flux
        %
        CN=min(F2(IJ),0.); % Mass fluxes, upwind from E
        CP=max(F2(IJ),0.); % Mass fluxes, upwind from P
        %
        FUUDS=CP*U(IJ)+CN*U(IJN);            % Explicit convective fluxes UDS, U
        FVUDS=CP*V(IJ)+CN*V(IJN);            % Explicit convective fluxes UDS, V
        FUCDS=F2(IJ)*(U(IJN)*FYN+U(IJ)*FYP); % Explicit convective fluxes CDS, U
        FVCDS=F2(IJ)*(V(IJN)*FYN+V(IJ)*FYP); % Explicit convective fluxes CDS, V       
        %
        AN(IJ )=+CN-D;                       % Coefficient for E
        AS(IJN)=-CP-D;                       % Coefficient for W
        %
        SU(IJ )=SU(IJ )+GDS*(FUUDS-FUCDS);   % Source term for P, U
        SU(IJN)=SU(IJN)-GDS*(FUUDS-FUCDS);   % Source term for N, U
        SV(IJ )=SV(IJ )+GDS*(FVUDS-FVCDS);   % Source term for P, V
        SV(IJN)=SV(IJN)-GDS*(FVUDS-FVCDS);   % Source term for N, V       
    end;
end;
        
%%% VOLUME INTEGRALS (SOURCE TERMS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=2:NIM,                                       % Loop I direction
    DX=X(I)-X(I-1);                                % Cell length
    for J=2:NJM,                                   % Loop J direction
        DY=Y(J)-Y(J-1);                            % Cell heigth
        VOL=DX*DY;                                 % Cell Volume
        IJ=LI(I)+J;                                % IJ index to P
        %
        PE=P(IJ+NJ)*FX(I  )+P(IJ   )*(1.-FX(I  )); % Pressure at 'e'
        PW=P(IJ   )*FX(I-1)+P(IJ-NJ)*(1.-FX(I-1)); % Pressure at 'w'
        PN=P(IJ+ 1)*FY(J  )+P(IJ   )*(1.-FY(J  )); % Pressure at 'n'
        PS=P(IJ   )*FY(J-1)+P(IJ- 1)*(1.-FY(J-1)); % Pressure at 's'
        DPX(IJ)=(PE-PW)/DX;                        % Pressure X-gradient
        DPY(IJ)=(PN-PS)/DY;                        % Pressure Y-gradient
        SU(IJ)=SU(IJ)-DPX(IJ)*VOL;                 % Pressure term - U
        SV(IJ)=SV(IJ)-DPY(IJ)*VOL;                 % Pressure term - V
        %
        if LTIME,                                  % Unsteady ?
            APT=DENSIT*VOL*DTR;                    % Coefficient
            SU(IJ)=SU(IJ)+APT*U0(IJ);              % Source term - U
            SV(IJ)=SV(IJ)+APT*V0(IJ);              % Source term - V
            AP(IJ)=AP(IJ)+APT;                     % AP Coef.
        end;
    end;
end;

%%% BOUNDARY CONDITIONS U,V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SOUTH BOUNDARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=2:NIM,                                  % Loop I direction
    IJ=LI(I)+2; IJB=IJ-1;                     % IJ,IJB to South Boundary (J=2)
    U(IJB)=0.;                                % Set South wall U velocity
    V(IJB)=0.;                                % Set South wall V velocity
    D=VISC*(X(I)-X(I-1))/(YC(2)-YC(1));       % Diffusion coef.
    AP(IJ)=AP(IJ)+D;                          % Coef for P
    SU(IJ)=SU(IJ)+D*U(IJB);                   % Sourec term - U
    SV(IJ)=SV(IJ)+D*V(IJB);                   % Sourec term - V   
end;

%%% NORTH BOUNDARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=2:NIM,
    IJ=LI(I)+NJM; IJB=IJ+1;
    U(IJB)=ULID;           
    V(IJB)=0.;
    D=VISC*(X(I)-X(I-1))/(YC(NJ)-YC(NJM));
    AP(IJ)=AP(IJ)+D;
    SU(IJ)=SU(IJ)+D*U(IJB);
    SV(IJ)=SV(IJ)+D*V(IJB);    
end;

%%% WEST BOUNDARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for J=2:NJM,
    IJ=LI(2)+J; IJB=IJ-NJ;
    U(IJB)=0.;
    V(IJB)=0.;
    D=VISC*(Y(J)-Y(J-1))/(XC(2)-XC(1));
    AP(IJ)=AP(IJ)+D;
    SU(IJ)=SU(IJ)+D*U(IJB);    
    SV(IJ)=SV(IJ)+D*V(IJB);
end;
     
%%% EAST BOUNDARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for J=2:NJM,
    IJ=LI(NIM)+J; IJB=IJ+NJ;
    U(IJB)=0.;
    V(IJB)=0.;
    D=VISC*(Y(J)-Y(J-1))/(XC(NI)-XC(NIM));
    AP(IJ)=AP(IJ)+D;
    SU(IJ)=SU(IJ)+D*U(IJB);    
    SV(IJ)=SV(IJ)+D*V(IJB);
end;

%%% SOLVE EQUATIONS FOR U AND V %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% UNDER RELAXATION, SOLVING FOR U-VELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=2:NIM,
    for IJ=LI(I)+2:LI(I)+NJM,
        APR(IJ)=AP(IJ);
        AP(IJ)=(APR(IJ)-AE(IJ)-AW(IJ)-AN(IJ)-AS(IJ))/URF(IU);
        SU(IJ)=SU(IJ)+(1.-URF(IU))*AP(IJ)*U(IJ);
    end;
end;

%%% SOLVE for U
U=SIPSOL(U,IU);

%%% UNDER RELAXATION, SOLVING FOR V-VELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=2:NIM,
    for IJ=LI(I)+2:LI(I)+NJM,
        AP(IJ)=(APR(IJ)-AE(IJ)-AW(IJ)-AN(IJ)-AS(IJ))/URF(IV);
        SU(IJ)=SV(IJ)+(1.-URF(IV))*AP(IJ)*V(IJ);
        APR(IJ)=1./AP(IJ);
    end;
end;

%%% SOLVE for V
V=SIPSOL(V,IV);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CALCP;

global LTIME LTEST
global IPR JPR IU IV IP
global NI NJ NIM NJM NIJ LI
global X Y XC YC FX FY
global DENSIT VISC GDS DT DTR SMALL ALFA ULID
global U V P PP F1 F2 DPX DPY U0 V0
global AE AW AS AN AP APR SU SV
global SOR URF NSW RESOR

%%% INITIALIZE COEFFICIENTS AND SOURCE TERMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AP,AE,AW,AN,AS,SU]=deal(zeros(NIJ,1));

%%% EAST CV FACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=2:NIM-1,              % Loop I direction
    DXPE=XC(I+1)-XC(I);     % Distance P->E
    FXE=FX(I);              % Interpolation factor
    FXP=1.-FXE;             % Interpolation factor
    %
    for J=2:NJM,            % Loop J direction
        IJ=LI(I)+J;         % Index IJ to P
        IJE=IJ+NJ;          % Index IJE to E
        S=Y(J)-Y(J-1);      % Cell face 'area'
        VOLE=DXPE*S;        % Volume between P and E
        D=DENSIT*S;         % Coefficient
        %
        DPXEL=0.5*(DPX(IJE)+DPX(IJ));   % Interpolated gradient in 'e'
        UEL=U(IJE)*FXE+U(IJ)*FXP;       % Interpolated U-velocity in 'e'
        APUE=APR(IJE)*FXE+APR(IJ)*FXP;  % Interpolated 1/AP coeff. in 'e'
        %
        DPXE=(P(IJE)-P(IJ))/DXPE;       % Gradient in 'e', compact aprox.
        UE=UEL-APUE*VOLE*(DPXE-DPXEL);  % Corrected U velocity in 'e'
        F1(IJ)=D*UE;                    % Mass flux through 'e' face
        %
        AE(IJ)=-D*APUE*S;               % Coefficient for E in P
        AW(IJE)=AE(IJ);                 % Coefficient for W in E
    end;
end;

%%% NORTH CV FACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for J=2:NJM-1,              % Loop J direction
    DYPN=YC(J+1)-YC(J);     % Distance P->N
    FYN=FY(J);              % Interpolation factor
    FYP=1.-FYN;             % Interpolation factor
    %
    for I=2:NIM,            % Loop I direction
        IJ=LI(I)+J;         % Index IJ to P
        IJN=IJ+1;           % Index IJN to N
        S=X(I)-X(I-1);      % Cell face 'area'
        VOLN=DYPN*S;        % Volume between P and N
        D=DENSIT*S;         % Coefficient
        %
        DPYNL=0.5*(DPY(IJN)+DPY(IJ));   % Interpolated gradient in 'n'
        VNL=V(IJN)*FYN+V(IJ)*FYP;       % Interpolated V-velocity in 'n'
        APVN=APR(IJN)*FYN+APR(IJ)*FYP;  % Interpolated 1/AP coeff. in 'n'
        %
        DPYN=(P(IJN)-P(IJ))/DYPN;       % Gradient in 'n', compact aprox.
        VN=VNL-APVN*VOLN*(DPYN-DPYNL);  % Corrected V velocity in 'n'
        F2(IJ)=D*VN;                    % Mass flux through 'n' face
        %
        AN(IJ)=-D*APVN*S;               % Coefficient for N in P
        AS(IJN)=AN(IJ);                 % Coefficient for S in N
    end;
end;

%%% SINCE ALL BOUNDARIES ARE ZERO MASS FLUX BOUNDARIES (WALLS),       %%%%%
%%% WE HAVE NEWMANN BOUNDARY CONDITIONS FOR PRESSURE CORRECTION       %%%%%
%%% (NULL GRADIENT). NO SPECIAL TREATMEN IS REQUIRED.                 %%%%%
%%% FOR THE CASE OF INLETS AND OUTLETS, MASS FLUXES AT THE BOUNDARIES %%%%%
%%% NEED TO BE COMPUTED HERE.                                         %%%%%

%%% SOURCE TERM AND COEFFICIENT FOR NODE P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUM=0.;                                           % Initialize SUM
for I=2:NIM,                                      % Loop I direction
    for IJ=LI(I)+2:LI(I)+NJM,                     % Loop J direction
        SU(IJ)=F1(IJ-NJ)-F1(IJ)+F2(IJ-1)-F2(IJ);  % Flux imbalance
        AP(IJ)=-(AE(IJ)+AW(IJ)+AN(IJ)+AS(IJ));    % Coefficient for P
        SUM=SUM+SU(IJ);                           % Global flux imbalance
        PP(IJ)=0.;                                % Initialize PP
    end;
end;

if LTEST, disp(['SUM=',num2str(SUM)]); end;       % If testing display SUM

%%% solve for PP
PP=SIPSOL(PP,IP);

%%% EXTRAPOLATE BOUNDARY VALUES FOR PP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOUTH AND NORTH BOUNDARIES
for I=2:NIM,
    IJ=LI(I)+1;
    PP(IJ)=PP(IJ+1)+(PP(IJ+1)-PP(IJ+2))*FY(2);
    IJ=LI(I)+NJ;
    PP(IJ)=PP(IJ-1)+(PP(IJ-1)-PP(IJ-2))*(1.-FY(NJM-1));
end;

%%% WEST AND EAST BOUNDARIES
NJ2=2*NJ;
for J=2:NJM,
    IJ=LI(1)+J;
    PP(IJ)=PP(IJ+NJ)+(PP(IJ+NJ)-PP(IJ+NJ2))*FX(2);
    IJ=LI(NI)+J;
    PP(IJ)=PP(IJ-NJ)+(PP(IJ-NJ)-PP(IJ-NJ2))*(1.-FX(NIM-1));
end;

%%% REFERENCE PRESSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IJPREF=LI(IPR)+JPR;
PP0=PP(IJPREF);

%%% CORRECT EAST MASS FLUXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=2:NIM-1
    for IJ=LI(I)+2:LI(I)+NJM,
        F1(IJ)=F1(IJ)+AE(IJ)*(PP(IJ+NJ)-PP(IJ));
    end;
end;

%%% CORRECT NORTH MASS FLUXES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I=2:NIM
    for IJ=LI(I)+2:LI(I)+NJM-1,
        F2(IJ)=F2(IJ)+AN(IJ)*(PP(IJ+1)-PP(IJ));
    end;
end;

%%% CORRECT PRESSURE AND VELOCITIES AT CELL CENTER %%%%%%%%%%%%%%%%%%%%%%%%
for I=2:NIM,
    DX=X(I)-X(I-1);        
    %
    for J=2:NJM,
        IJ=LI(I)+J;
        DY=Y(J)-Y(J-1);
        %
        PPE=PP(IJ+NJ)*FX(I  )+PP(IJ   )*(1.-FX(I  )); % Pressure at 'e'
        PPW=PP(IJ   )*FX(I-1)+PP(IJ-NJ)*(1.-FX(I-1)); % Pressure at 'w'
        PPN=PP(IJ+1 )*FY(J  )+PP(IJ   )*(1.-FY(J  )); % Pressure at 'n'
        PPS=PP(IJ   )*FY(J-1)+PP(IJ-1 )*(1.-FY(J-1)); % Pressure at 's'
        %
        U(IJ)=U(IJ)-(PPE-PPW)*DY*APR(IJ);  % U-velocity corrected
        V(IJ)=V(IJ)-(PPN-PPS)*DX*APR(IJ);  % V-velocity corrected
        P(IJ)=P(IJ)+URF(IP)*(PP(IJ)-PP0);  % Pressure corrected
    end;
end;

%%% EXTRAPOLATE BOUNDARY VALUES FOR P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOUTH AND NORTH BOUNDARIES
for I=2:NIM,
    IJ=LI(I)+1;
    P(IJ)=P(IJ+1)+(P(IJ+1)-P(IJ+2))*FY(2);
    IJ=LI(I)+NJ;
    P(IJ)=P(IJ-1)+(P(IJ-1)-P(IJ-2))*(1.-FY(NJM-1));
end;

%%% WEST AND EAST BOUNDARIES
NJ2=2*NJ;
for J=2:NJM,
    IJ=LI(1)+J;
    P(IJ)=P(IJ+NJ)+(P(IJ+NJ)-P(IJ+NJ2))*FX(2);
    IJ=LI(NI)+J;
    P(IJ)=P(IJ-NJ)+(P(IJ-NJ)-P(IJ-NJ2))*(1.-FX(NIM-1));
end;

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FI=SIPSOL(FI,IFI);

global LTEST
global NI NJ NIM NJM NIJ LI
global SMALL ALFA
global AE AW AS AN AP SU
global SOR NSW RESOR

%%% INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[UE,UN,RES,LW,LS,LPR]=deal(zeros(NIJ,1));

%%% COEFFICIENTS OF UPPER AND LOWER TIRANGULAR MATRICES %%%%%%%%%%%%%%%%%%%
for I=2:NIM,
    for IJ=LI(I)+2:LI(I)+NJM,
        LW(IJ)=AW(IJ)/(1.+ALFA*UN(IJ-NJ));
        LS(IJ)=AS(IJ)/(1.+ALFA*UE(IJ-1));
        P1=ALFA*LW(IJ)*UN(IJ-NJ);
        P2=ALFA*LS(IJ)*UE(IJ-1);
        LPR(IJ)=1./(AP(IJ)+P1+P2-LW(IJ)*UE(IJ-NJ)-LS(IJ)*UN(IJ-1));
        UN(IJ)=(AN(IJ)-P1)*LPR(IJ);
        UE(IJ)=(AE(IJ)-P2)*LPR(IJ);
    end;
end;

%%% INNER ITERATIONS LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for L=1:NSW(IFI),
    RESL=0.;
    %
    %%% CALCULATE RESIDUAL AND OVERWRITE IT BY INTERMEDIATE VECTOR%%%%%%%%%
    for I=2:NIM,
        for IJ=LI(I)+2:LI(I)+NJM,
            RES(IJ)=SU(IJ)-AN(IJ)*FI(IJ+1)-AS(IJ)*FI(IJ-1)-AE(IJ)*FI(IJ+NJ)...
                   -AW(IJ)*FI(IJ-NJ)-AP(IJ)*FI(IJ);
            RESL=RESL+abs(RES(IJ));
            RES(IJ)=(RES(IJ)-LS(IJ)*RES(IJ-1)-LW(IJ)*RES(IJ-NJ))*LPR(IJ);
        end;
    end;
    %
    %%% STORE INITIAL RESIDUAL SUM FOR CHECKING CONVERGENCE OF OUTER IT. %%
    if L==1, RESOR(IFI)=RESL; end;
    RSM=RESL/(RESOR(IFI)+SMALL);
    %
    %%% BACK SUBSTITUTION AND CORRECTION
    for I=NIM:-1:2,
        for IJ=LI(I)+NJM:-1:LI(I)+2,
            RES(IJ)=RES(IJ)-UN(IJ)*RES(IJ+1)-UE(IJ)*RES(IJ+NJ);
            FI(IJ)=FI(IJ)+RES(IJ);
        end;
    end;
    %
    %%% CHECK CONVERGENCE
    if LTEST, disp([num2str(L),' INNER ITER., RESL=',num2str(RESL)]); end;
    if RSM<SOR(IFI), return; end;
end;
%
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%