function varargout = vonNeumRichtVisc(FF, GR, BC, varargin)

alpha = 0.01;%005;
% power = 2;
% visc = zeros(size(FF, 1), size(FF,2));
epsMin = 0.005;
% visc = zeros(size(FF));
% indRho = reshape(strcmp(BC.N.varName, '\rho'),1,1,size(FF,3));
% indV1 = reshape(strcmp(BC.N.varType, 'v1'),1,1,size(FF,3));
% indV2 = reshape(strcmp(BC.N.varType, 'v2'),1,1,size(FF,3));
% bcVals = bcCalc(GR,BC,FF,DIR,find(indRho));
% indV1 = ones(1,1,size(FF,3));
% indV2 = zeros(size(indV1));
DIRnames = {'N', 'S', 'E', 'W'};

if ~isempty(varargin)
    varargout = calcVisc(varargin{1});
else
    for i = 1:4
        varargout{i} = calcVisc(DIRnames{i});
    end
end


function visc = calcVisc(DIR)
    bcVals = bcCalc(GR,BC,FF,DIR);

    switch DIR
        case 'W' 
            rhoSens = 1;%0.5.*(FF(:,:,indRho) + [bcVals(:,:,indRho), FF(:,1:end-1, indRho)]);
    %         dudx = 
            visc = alpha .*GR.dx .* rhoSens.* abs(diff([bcVals, FF], 1, 2)).^2 ./ (0.5.*abs(FF + [bcVals, FF(:,1:end-1,:)])) ;

        case 'E'
            rhoSens = 1;%0.5.*(FF(:,:,indRho) + [FF(:,2:end, indRho), bcVals(:,:,indRho)]);
            visc = alpha .*GR.dx.* rhoSens.* abs(diff([FF, bcVals], 1, 2)).^2 ./ (0.5.*abs(FF + [FF(:,2:end, :), bcVals])) ;

        case 'N'
            rhoSens = 1;%0.5.*(FF(:,:,indRho) + [FF(2:end,:,indRho); bcVals(:,:,indRho)]);
            visc = alpha .*GR.dy.* rhoSens.* abs(diff([FF; bcVals], 1, 1)).^2 ./ (0.5.*abs(FF + [FF(2:end,:,:); bcVals])) ;

        case 'S'
            rhoSens = 1;%0.5.*(FF(:,:,indRho) + [bcVals(:,:,indRho); FF(1:end-1,:, indRho)]);
            visc = alpha .*GR.dy.* rhoSens.* abs(diff([bcVals; FF], 1, 1)).^2 ./ (0.5.*abs(FF + [bcVals; FF(1:end-1,:,:)])) ;

    end
    testInd = isinf(visc) + isnan(visc);
    visc(logical(testInd)) = 0;
    visc = visc + epsMin;
end

% visc = repmat(visc, 1,1,size(FF,3));
% if max(visc(:)) > 0
%     visc = repmat(visc./max(visc(:)), 1,1,size(FF,3));
% else
%     visc = repmat(visc, 1,1,size(FF,3));
% end

end