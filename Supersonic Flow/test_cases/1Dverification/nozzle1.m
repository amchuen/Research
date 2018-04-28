clc;
clear;
close all;

%% Define Points

x=[0,.5,1,1.5,2,2.5,3,3.5,4,4.5,4.6,4.7,4.8,4.9,5,5.1,5.2,...
    5.3,5.4,5.5,5.5957,5.6823,5.7812,5.9096,6.0159,...
    6.1438,6.2892,6.4219,6.5755,6.6636,6.8324,6.9493,7.0804,...
    7.3283,7.5008,7.5807,7.7265,7.8414,8.0561,8.3488,8.5615,...
    8.9019,9.2759,9.6865,10.1373,10.6322,11.1753,11.7716,...
    12.4262,13.1447,13.9339,14.8006,15.7529,16.7996,17.9504,...
    19.2166,20.0316];

y=[2.636,2.4097,2.1848,1.9618,1.7414,1.5248,1.3138,1.1117,...
    0.9241,0.7616,0.7336,0.7075,0.6835,0.6619,0.6428,...
    0.6265,0.6131,0.6029,0.5960,0.5924,0.5962,0.6033,...
    0.6122,0.6289,0.6478,0.6747,0.7085,0.7434,0.7836,0.8082,...
    0.8574,0.8924,0.9330,1.0104,1.0645,1.0897,1.1357,...
    1.1719,1.2395,1.3294,1.3928,1.4897,1.5893,1.6911,1.7946,...
    1.8992,2.0042,2.1088,2.2119,2.3123,2.4086,2.4992,...
    2.5819,2.6545,2.7142,2.7577,2.7759];

% Plot
plot(x,y, 'o')
axis('equal')
axis([0,25,-3,3])
hold on
plot(x,-y, 'o')
a=2*y;
% scatter(x,y)
% scatter(x,-y)
xlabel('x [inches]')
ylabel('y [inches]')

%% Generate Spline

Area = 2.*y;
curve = spline(x./min(Area), Area./min(Area));

xEnd = x(end)./min(Area);
xx = linspace(0,xEnd);
yy = ppval(curve,xx);

hold on; plot(xx.*min(Area), 0.5.*yy.*min(Area));

func = @(xval) ppval(curve, xval);

options = optimset('TolX', 1e-10);
xThroat = fminbnd(func, 0, xEnd, options);

save('mcCabe_nozzle.mat', 'curve', 'xThroat', 'xEnd');