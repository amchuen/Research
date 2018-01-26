%%% Run pcolm
%[X,Y,XC,YC,FX,FY,U,V,P,F1,F2]=pcolm;

%%% Get size
NI=length(X); NJ=length(Y); sz=[NI,NJ];

%%% Quiver velocity field
figure;
quiver(XC,YC,reshape(U,sz),reshape(V,sz));
axis image;

%%% Streamlines
figure;
h=streamline(XC',YC,reshape(U,sz),reshape(V,sz),rand(50,1),rand(50,1),[0.1,1000]);
axis image;
axis([0,1,0,1]);