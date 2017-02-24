clear;
close all;
clc;

%% Setup

% Flow Conditions
M0 = 1.5;
Beta = sqrt(abs(M0^2 - 1));

% Domain
dx = 0.01;
dy_bar = dx;
xvals = -10*dx:dx:(2+dx);
ybar = 0:dy_bar:4;
yvals = ybar./Beta;
dy = yvals(2) - yvals(1);

[XX, YY] = meshgrid(xvals, yvals);
% [~, YYbar] = meshgrid(xvals, ybar);

%% Body Geometry
% xB = 0:dx:1;
% yB = [(0:dx:0.5).*0.25, (- 0.25*0.5)/(1 - 0.5) .* (((0.5+dx):dx:1) - 0.5)];
yB = zeros(size(xvals));

% Diamond Body
yB(xvals>0 & xvals<=0.5) = xvals((xvals>0 & xvals<=0.5)).*0.25;
yB(xvals>(0.5) & xvals<=1) = (-0.25*0.5)/(1-0.5) .* (xvals(xvals>(0.5) & xvals<=1)-0.5) + 0.25*0.5;

% Parabolic Body
% tau = 0.25;
% yB(xvals>=0 & xvals <=1) = 2*tau * xvals(xvals>=0 & xvals <=1) .*(1 - xvals(xvals>=0 & xvals <=1));

%% Numerical Solution
% Potential Field
PHI = sspf_explicit(XX, YY, yB, M0);

% Velocity Field
UVALS = ones(length(yvals),length(xvals));
VVALS = zeros(size(UVALS));
for ii = 1:length(xvals)
    if ii == 1
        UVALS(:,ii) = ones(size(UVALS(:,ii)));
    elseif ii == length(xvals)
        UVALS(:,ii) = ones(size(UVALS(:,ii)));
    else
        UVALS(:,ii) = (PHI(:,ii+1) - PHI(:,ii-1))./(2*dx);
    end
end

% for jj = 2:length(yvals) - 1
%    VVALS(jj,:) = (PHI(jj+1,1:end-1) - PHI(jj-1,1:end-1))./(2*dy); 
% end
% 
% UVALS = round(UVALS, 5); % round results to remove small numerical noise
% VVALS = round(VVALS, 5);

figure();
contourf(XX(:,1:end-1),YY(:,1:end-1),PHI(:,1:end-1), 20);

figure();
contourf(XX,YY,UVALS, 20);
colorbar('eastoutside');

figure();
plot(xvals,UVALS(1,:));

% 
% figure();
% contourf(XX(:,1:end-1),YY(:,1:end-1),hypot(UVALS,VVALS), 20);
% colorbar('eastoutside');