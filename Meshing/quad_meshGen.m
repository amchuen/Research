clc;
clear
close all;

%% Initialize Computational Grid and Physical Grid

% Nozzle Function
nozzCoeffs = [2 -2 1]';

% Computational Grid
nZeta = 11; % dZ = 1
nEta = 21; % dE = 1
[ZZ, EE] = meshgrid(1:nZeta, 1:nEta);

% Physical Grid -> boundary conditions
xRange = [0.5, 1];
xS = equalCurveDist(nozzCoeffs, xRange, 2*(nZeta-1)+1);
xN = xS;
xW = linspace(xS(1), xN(1), 2*(nEta-1)+1);%, 1);
xE = ones(2*(nEta-1)+1, 1);

yN = [0.5.*(1 + (2.*xN(xN <1)-1).^2), ones(size(xN(xN >=1)))];
yS = -yN;
yW = linspace(yS(1), yN(1), 2*(nEta-1)+1);
yE = linspace(yS(end), yN(end), 2*(nEta-1)+1);

% Generate Primary Grid
UU_1 = cat(3, ZZ, EE); % [X, Y]
normU = zeros(length(yS),2);
normL = normU;
coeffMat_1 = zeros(4,size(UU_1,2));
for i = 1:size(xN,2)
%     ind = 2*i-1;
%     UU(:,i,1:2) = cat(3, linspace(xS(ind), xN(ind), nEta), linspace(yS(ind), yN(ind), nEta));
%     if i == 1
%         normU(i,:) = [-(-1.5*yN(i)+2*yN(i+1)-0.5*yN(i+2)),(-1.5*xN(i)+2*xN(i+1)-0.5*xN(i+2))]./sqrt((-1.5*yN(i)+2*yN(i+1)-0.5*yN(i+2))^2+(-1.5*xN(i)+2*xN(i+1)-0.5*xN(i+2))^2);
%         normL(i,:) = [-(-1.5*yS(i)+2*yS(i+1)-0.5*yS(i+2)),(-1.5*xS(i)+2*xS(i+1)-0.5*xS(i+2))]./sqrt((-1.5*yS(i)+2*yS(i+1)-0.5*yS(i+2))^2+(-1.5*xS(i)+2*xS(i+1)-0.5*xS(i+2))^2);
%     elseif i == size(xN,2) %size(UU,2)
%         normU(i,:) = [-(1.5*yN(i)-2*yN(i-1)+0.5*yN(i-2)),(1.5*xN(i)-2*xN(i-1)+0.5*xN(i-2))]./sqrt((1.5*yN(i)-2*yN(i-1)+0.5*yN(i-2)).^2 + (1.5*xN(i)-2*xN(i-1)+0.5*xN(i-2)).^2);
%         normL(i,:) = [-(1.5*yS(i)-2*yS(i-1)+0.5*yS(i-2)),(1.5*xS(i)-2*xS(i-1)+0.5*xS(i-2))]./sqrt((1.5*yS(i)-2*yS(i-1)+0.5*yS(i-2)).^2 + (1.5*xS(i)-2*xS(i-1)+0.5*xS(i-2)).^2);
%     else
%         normU(i,:) = [-(yN(i+1)-yN(i-1)),(xN(i+1)-xN(i-1))]./sqrt((yN(i+1)-yN(i-1))^2+(xN(i+1)-xN(i-1))^2);
%         normL(i,:) = [-(yS(i+1)-yS(i-1)),(xS(i+1)-xS(i-1))]./sqrt((yS(i+1)-yS(i-1))^2+(xS(i+1)-xS(i-1))^2);
%     end
    normU(i,:) = [-2*(2*xN(i)-1), 1]./sqrt((-2*(2*xN(i)-1)).^2 + 1);
    normL(i,:) = [2*(2*xN(i)-1), 1]./sqrt((2*(2*xN(i)-1)).^2 + 1);
    
    if mod(i,2) ~= 0 % update primary grid
        ind = 0.5*(i + 1);
        Amat = [ yS(i)^3, yS(i)^2, yS(i), 1; 3*yS(i)^2, 2*(yS(i)), 1, 0;...
                yN(i)^3, yN(i)^2, yN(i), 1; 3*yN(i)^2, 2*(yN(i)), 1, 0];
        bVec = [xS(i), normL(i,1)/normL(i,2), xN(i), normU(i,1)/normU(i,2)]';
        coeffMat_1(:,ind) = Amat\bVec;
        UU_1(:,ind,2) = equalCurveDist(coeffMat_1(:,ind), [yS(i), yN(i)], length(UU_1(:,ind,2)));
        UU_1(:,ind,1) = polyval(coeffMat_1(:,ind), UU_1(:,ind,2));
    end

end

% build secondary grid
UU_2 = zeros(size(UU_1));
UU_2(end+1, end+1, :) = 0;
UU_2(:, 1,1) = xW([1,2:2:length(xW)-1,length(xW)]);
UU_2(:, 1,2) = yW([1,2:2:length(xW)-1,length(xW)]);
UU_2(end, :,1) = xN([1,2:2:length(xN)-1,length(xN)]);
UU_2(end, :,2) = yN([1,2:2:length(yN)-1,length(yN)]);
UU_2(1, :,1) = xS([1,2:2:length(xS)-1,length(xS)]);
UU_2(1, :,2) = yS([1,2:2:length(yS)-1,length(yS)]);
coeffMat_2 = zeros(4,size(UU_2,2));
for i = [1,2:2:length(xN)-1,length(xN)] %2:2:(size(xN,2)-1)
    if i == 1
        Amat = [ yS(i)^3, yS(i)^2, yS(i), 1; 3*yS(i)^2, 2*(yS(i)), 1, 0;...
                yN(i)^3, yN(i)^2, yN(i), 1; 3*yN(i)^2, 2*(yN(i)), 1, 0];
        bVec = [xS(i), normL(i,1)/normL(i,2), xN(i), normU(i,1)/normU(i,2)]';
        coeffMat_2(:,i) = Amat\bVec;
%         UU_2(:,i,2) = equalCurveDist(coeffMat_2(:,i), [yS(i), yN(i)], length(UU_2(:,i,2)));
%         UU_2(:,i,1) = polyval(coeffMat_2(:,i), UU_2(:,i,2));
    elseif i == length(xN)
        y_S = 0.5*(UU_1(1,end,2) + UU_1(2,end,2));
        y_N = 0.5*(UU_1(end,end,2) + UU_1(end-1,end,2));
        Amat = [ yS(i)^3, yS(i)^2, yS(i), 1; 3*yS(i)^2, 2*(yS(i)), 1, 0;...
                yN(i)^3, yN(i)^2, yN(i), 1; 3*yN(i)^2, 2*(yN(i)), 1, 0];
        bVec = [xS(i), normL(i,1)/normL(i,2), xN(i), normU(i,1)/normU(i,2)]';
        coeffMat_2(:,end) = Amat\bVec;
        UU_2(2:end-1,end,2) = equalCurveDist(coeffMat_2(:,end), [y_S, y_N], length(UU_2(2:end-1,end,2)));
        UU_2(2:end-1,end,1) = polyval(coeffMat_2(:,end), UU_2(2:end-1,end,2));
    else
        ind = 0.5*(i);
        y_S = 0.25*(UU_1(1,ind,2) + UU_1(1,ind+1,2) + UU_1(2,ind,2) + UU_1(2,ind+1,2));
        y_N = 0.25*(UU_1(end,ind,2) + UU_1(end,ind+1,2) + UU_1(end-1,ind,2) + UU_1(end-1,ind+1,2));
        Amat = [ yS(i)^3, yS(i)^2, yS(i), 1; 3*yS(i)^2, 2*(yS(i)), 1, 0;...
                yN(i)^3, yN(i)^2, yN(i), 1; 3*yN(i)^2, 2*(yN(i)), 1, 0];
        bVec = [xS(i), normL(i,1)/normL(i,2), xN(i), normU(i,1)/normU(i,2)]';
        coeffMat_2(:,ind+1) = Amat\bVec;
        UU_2(2:end-1,ind+1,2) = equalCurveDist(coeffMat_2(:,ind+1), [y_S, y_N], length(UU_2(2:end-1,ind+1,2)));
        UU_2(2:end-1,ind+1,1) = polyval(coeffMat_2(:,ind+1), UU_2(2:end-1,ind+1,2));
    end
    
end

%% Fit a quadratic Across X

% Refit Primary Grid
normB_1 = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat_1(1:end-1,1),UU_1(:,1,2)).^2), polyval((3:-1:1)'.*coeffMat_1(1:end-1,1),UU_1(:,1,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat_1(1:end-1,1),UU_1(:,1,2)).^2)];
normMat_1 = zeros(size(UU_1,1), size(UU_1,2));
for i = 2:(size(UU_1,1)-1) % march in y-direction
    
    % Grab Normals at Starting Point
    y0p = normB_1(i,2)./normB_1(i,1);
    
    for ii = 2:size(UU_1,2) % march in x-direction
        func = @(yVec)  [   yVec(2) + polyval((3:-1:1)'.*coeffMat_1(1:end-1,ii), yVec(1));...
                            0.5.*(yVec(2)+y0p).*(polyval(coeffMat_1(:,ii),yVec(1))-UU_1(i,ii-1,1))+UU_1(i,ii-1,2)-yVec(1)];
                        
        jacob = @(yVec) [   polyval((2:-1:1)'.*(3:-1:2)'.*coeffMat_1(1:end-2,ii), yVec(1)), 1;...
                            0.5.*(yVec(2)+y0p).*polyval(coeffMat_1(1:end-1,ii),yVec(1))-1, 0.5.*(polyval(coeffMat_1(:,ii),yVec(1))-UU_1(i,ii-1,1))];
        [yNew, res] = newtonSys(func, jacob, [UU_1(i,ii,2); -polyval((3:-1:1)'.*coeffMat_1(1:end-1,ii),UU_1(i,1,2))], 1e-7, 1e-7, 1e5, 1);
        UU_1(i,ii,1) = polyval(coeffMat_1(:,ii),yNew(1)); UU_1(i,ii,2) = yNew(1);
        y0p = yNew(2);
        tVec = [1, y0p];
        nVec = [-1, polyval((3:-1:1)'.*coeffMat_1(1:end-1,ii),UU_1(i,ii,2))];
        normMat_1(i,ii) = 1-abs((tVec./norm(tVec))*(nVec./norm(nVec))');
    end
    
end

% Refit Secondary Grid
normB_2 = [-1./sqrt(1+polyval((3:-1:1)'.*coeffMat_2(1:end-1,1),UU_1(:,1,2)).^2), polyval((3:-1:1)'.*coeffMat_2(1:end-1,1),UU_1(:,1,2))./sqrt(1+polyval((3:-1:1)'.*coeffMat_2(1:end-1,1),UU_1(:,1,2)).^2)];
normMat_2 = zeros(size(UU_2,1), size(UU_2,2));
for i = 2:(size(UU_2,1)-1) % march in y-direction
    
    % Grab Normals at Starting Point
    y0p = normB_1(i,2)./normB_1(i,1);
    
    for ii = 2:size(UU_2,2) % march in x-direction
        func = @(yVec)  [   yVec(2) + polyval((3:-1:1)'.*coeffMat_2(1:end-1,ii), yVec(1));...
                            0.5.*(yVec(2)+y0p).*(polyval(coeffMat_2(:,ii),yVec(1))-UU_2(i,ii-1,1))+UU_2(i,ii-1,2)-yVec(1)];
                        
        jacob = @(yVec) [   polyval((2:-1:1)'.*(3:-1:2)'.*coeffMat_2(1:end-2,ii), yVec(1)), 1;...
                            0.5.*(yVec(2)+y0p).*polyval(coeffMat_2(1:end-1,ii),yVec(1))-1, 0.5.*(polyval(coeffMat_2(:,ii),yVec(1))-UU_2(i,ii-1,1))];
        [yNew, res] = newtonSys(func, jacob, [UU_2(i,ii,2); -polyval((3:-1:1)'.*coeffMat_2(1:end-1,ii),UU_2(i,1,2))], 1e-7, 1e-7, 1e5, 1);
        UU_2(i,ii,1) = polyval(coeffMat_2(:,ii),yNew(1)); UU_2(i,ii,2) = yNew(1);
        y0p = yNew(2);
        tVec = [1, y0p];
        nVec = [-1, polyval((3:-1:1)'.*coeffMat_2(1:end-1,ii),UU_2(i,ii,2))];
        normMat_2(i,ii) = 1-abs((tVec./norm(tVec))*(nVec./norm(nVec))');
    end
    
end

% close all;
figure(10); 
plot(UU_1(:,:,1), UU_1(:,:,2), 'b-', UU_1(:,:,1)', UU_1(:,:,2)', 'b-'); axis equal;hold on;
plot(UU_2(:,:,1), UU_2(:,:,2), 'r--', UU_2(:,:,1)', UU_2(:,:,2)', 'r--'); axis equal;
saveas(gcf, 'quadratic_spline_mesh', 'pdf');

%% Find Cell Volumes
UU_2SW = UU_2(1:end-1,1:end-1,:);
UU_2SE = UU_2(1:end-1,2:end,:);
UU_2NE = UU_2(2:end,2:end,:);
UU_2NW = UU_2(2:end,1:end-1,:);

% Calculate initial volumes first
UU_vol = 0.5.*abs(  UU_2SW(:,:,1).*UU_2SE(:,:,2) + ... (1,1)*(2,1)
                    UU_2SE(:,:,1).*UU_2NE(:,:,2) + ... (2,1)*(2,2)
                    UU_2NE(:,:,1).*UU_2NW(:,:,2) + ... (2,2)*(1,2)
                    UU_2NW(:,:,1).*UU_2SW(:,:,2) + ... (1,2)*(1,1)
                    -(UU_2SW(:,:,2).*UU_2SE(:,:,1) + ... (1,1)*(2,1)
                    UU_2SE(:,:,2).*UU_2NE(:,:,1) + ... (2,1)*(2,2)
                    UU_2NE(:,:,2).*UU_2NW(:,:,1) + ... (2,2)*(1,2)
                    UU_2NW(:,:,2).*UU_2SW(:,:,1)) ... (1,2)*(1,1)
                    );
                
% Correct Boundaries
UU_vol(2:end-1,end) = UU_vol(2:end-1,end)   - UU_2(2:end-2,end,1).*UU_2(3:end-1,end,2) + UU_2(2:end-2,end,2).*UU_2(3:end-1,end,1)... remove contribution of red edge at exit, need to add interpolations between red and blue edges!
                                            + UU_2(2:end-2,end,1).*UU_1(2:end-1,end,2) + UU_1(2:end-1,end,1).*UU_2(3:end-1,end,2)... add contributions of x in primary grid
                                            - UU_2(2:end-2,end,2).*UU_1(2:end-1,end,1) - UU_1(2:end-1,end,2).*UU_2(3:end-1,end,1); % subtract contributions of y in primary grid 
    
%% Calculate Grid Normals

normN = cat(3, -diff(UU_2(2:end,:,2),1,2), diff(UU_2(2:end,:,1),1,2));
normS = cat(3, -diff(UU_2(1:end-1,:,2),1,2), diff(UU_2(1:end-1,:,1),1,2));
normW = cat(3, -diff(UU_2(:,1:end-1,2),1,1), diff(UU_2(:,1:end-1,1),1,1));
normE = cat(3, -diff(UU_2(:,2:end,2),1,1), diff(UU_2(:,2:end,1),1,1));

%% Calculate Weights

UU_N = [UU_1(2:end,:,:);UU_1(end,:,:)];
UU_S = [UU_1(1,:,:);UU_1(1:end-1,:,:)];
UU_W = [UU_1(:,1,:),UU_1(:,1:end-1,:)];
UU_E = [UU_1(:,2:end,:),UU_1(:,end,:)];

midpt_N = cat(3,    ((UU_1(:,:,1).*UU_N(:,:,2) - UU_1(:,:,2).*UU_N(:,:,1)).*(UU_2NW(:,:,1)-UU_2NE(:,:,1))-(UU_1(:,:,1) - UU_N(:,:,1)).*(UU_2NW(:,:,1).*UU_2NE(:,:,2)-UU_2NW(:,:,2).*UU_2NE(:,:,1)))./...
                    ((UU_1(:,:,1)-UU_N(:,:,1)).*(UU_2NW(:,:,2) - UU_2NE(:,:,2)) - (UU_1(:,:,2)-UU_N(:,:,2)).*(UU_2NW(:,:,1)-UU_2NE(:,:,1))),...
                    ((UU_1(:,:,1).*UU_N(:,:,2) - UU_1(:,:,2).*UU_N(:,:,1)).*(UU_2NW(:,:,2)-UU_2NE(:,:,2))-(UU_1(:,:,2) - UU_N(:,:,2)).*(UU_2NW(:,:,1).*UU_2NE(:,:,2)-UU_2NW(:,:,2).*UU_2NE(:,:,1)))./...
                    ((UU_1(:,:,1)-UU_N(:,:,1)).*(UU_2NW(:,:,2) - UU_2NE(:,:,2)) - (UU_1(:,:,2)-UU_N(:,:,2)).*(UU_2NW(:,:,1)-UU_2NE(:,:,1)))...
                    );
midpt_N(end,:,:) = UU_1(end,:,:);

wt_N = cat(3, sqrt(sum((UU_1 - midpt_N).^2,3)), sqrt(sum((UU_N - midpt_N).^2,3)));
wt_N(1:end-1,:,:) = wt_N(1:end-1,:,:)./sqrt(sum((UU_1(1:end-1,:,:)-UU_N(1:end-1,:,:)).^2,3));

midpt_S = cat(3,    ((UU_1(:,:,1).*UU_S(:,:,2) - UU_1(:,:,2).*UU_S(:,:,1)).*(UU_2SW(:,:,1)-UU_2SE(:,:,1))-(UU_1(:,:,1) - UU_S(:,:,1)).*(UU_2SW(:,:,1).*UU_2SE(:,:,2)-UU_2SW(:,:,2).*UU_2SE(:,:,1)))./...
                    ((UU_1(:,:,1)-UU_S(:,:,1)).*(UU_2SW(:,:,2) - UU_2SE(:,:,2)) - (UU_1(:,:,2)-UU_S(:,:,2)).*(UU_2SW(:,:,1)-UU_2SE(:,:,1))),...
                    ((UU_1(:,:,1).*UU_S(:,:,2) - UU_1(:,:,2).*UU_S(:,:,1)).*(UU_2SW(:,:,2)-UU_2SE(:,:,2))-(UU_1(:,:,2) - UU_S(:,:,2)).*(UU_2SW(:,:,1).*UU_2SE(:,:,2)-UU_2SW(:,:,2).*UU_2SE(:,:,1)))./...
                    ((UU_1(:,:,1)-UU_S(:,:,1)).*(UU_2SW(:,:,2) - UU_2SE(:,:,2)) - (UU_1(:,:,2)-UU_S(:,:,2)).*(UU_2SW(:,:,1)-UU_2SE(:,:,1)))...
                    );
midpt_S(1,:,:) = UU_1(1,:,:);
wt_S = cat(3, sqrt(sum((UU_1 - midpt_S).^2,3)), sqrt(sum((UU_S - midpt_S).^2,3)));
wt_S(2:end,:,:) = wt_S(2:end,:,:)./sqrt(sum((UU_1(2:end,:,:)-UU_S(2:end,:,:)).^2,3));


midpt_E = cat(3,    ((UU_1(:,:,1).*UU_E(:,:,2) - UU_1(:,:,2).*UU_E(:,:,1)).*(UU_2SE(:,:,1)-UU_2NE(:,:,1))-(UU_1(:,:,1) - UU_E(:,:,1)).*(UU_2SE(:,:,1).*UU_2NE(:,:,2)-UU_2SE(:,:,2).*UU_2NE(:,:,1)))./...
                    ((UU_1(:,:,1)-UU_E(:,:,1)).*(UU_2SE(:,:,2) - UU_2NE(:,:,2)) - (UU_1(:,:,2)-UU_E(:,:,2)).*(UU_2SE(:,:,1)-UU_2NE(:,:,1))),...
                    ((UU_1(:,:,1).*UU_E(:,:,2) - UU_1(:,:,2).*UU_E(:,:,1)).*(UU_2SE(:,:,2)-UU_2NE(:,:,2))-(UU_1(:,:,2) - UU_E(:,:,2)).*(UU_2SE(:,:,1).*UU_2NE(:,:,2)-UU_2SE(:,:,2).*UU_2NE(:,:,1)))./...
                    ((UU_1(:,:,1)-UU_E(:,:,1)).*(UU_2SE(:,:,2) - UU_2NE(:,:,2)) - (UU_1(:,:,2)-UU_E(:,:,2)).*(UU_2SE(:,:,1)-UU_2NE(:,:,1)))...
                    );
midpt_E(:,end,:) = UU_1(:,end,:);
wt_E = cat(3, sqrt(sum((UU_1 - midpt_E).^2,3)), sqrt(sum((UU_E - midpt_E).^2,3)));
wt_E(:,1:end-1,:) = wt_E(:,1:end-1,:)./sqrt(sum((UU_1(:,1:end-1,:)-UU_E(:,1:end-1,:)).^2,3));


midpt_W = cat(3,    ((UU_1(:,:,1).*UU_W(:,:,2) - UU_1(:,:,2).*UU_W(:,:,1)).*(UU_2NW(:,:,1)-UU_2SW(:,:,1))-(UU_1(:,:,1) - UU_W(:,:,1)).*(UU_2NW(:,:,1).*UU_2SW(:,:,2)-UU_2NW(:,:,2).*UU_2SW(:,:,1)))./...
                    ((UU_1(:,:,1)-UU_W(:,:,1)).*(UU_2NW(:,:,2) - UU_2SW(:,:,2)) - (UU_1(:,:,2)-UU_W(:,:,2)).*(UU_2NW(:,:,1)-UU_2SW(:,:,1))),...
                    ((UU_1(:,:,1).*UU_W(:,:,2) - UU_1(:,:,2).*UU_W(:,:,1)).*(UU_2NW(:,:,2)-UU_2SW(:,:,2))-(UU_1(:,:,2) - UU_W(:,:,2)).*(UU_2NW(:,:,1).*UU_2SW(:,:,2)-UU_2NW(:,:,2).*UU_2SW(:,:,1)))./...
                    ((UU_1(:,:,1)-UU_W(:,:,1)).*(UU_2NW(:,:,2) - UU_2SW(:,:,2)) - (UU_1(:,:,2)-UU_W(:,:,2)).*(UU_2NW(:,:,1)-UU_2SW(:,:,1)))...
                    );
midpt_W(:,1,:) = UU_1(:,1,:);
wt_W = cat(3, sqrt(sum((UU_1 - midpt_W).^2,3)), sqrt(sum((UU_W - midpt_W).^2,3)));
wt_W(:,2:end,:) = wt_W(:,2:end,:)./sqrt(sum((UU_1(:,2:end,:)-UU_W(:,2:end,:)).^2,3));

%% Save info to .mat file

save('quad_nozzle.mat', 'UU_1', 'UU_2', 'normN', 'normS', 'normE', 'normW',...
                        'midpt_N', 'midpt_S', 'midpt_E', 'midpt_W',...
                        'wt_N', 'wt_S', 'wt_E', 'wt_W',...
                        'UU_vol');