function varargout = vonNeumRichtVisc(FF, GR, BC, FL, varargin)

alpha = 0.25;
% power = 2;
epsMin = 0.25*GR.dx*GR.dy;
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
%             rhoSens = 1;%0.5.*(FF(:,:,indRho) + [bcVals(:,:,indRho), FF(:,1:end-1, indRho)]);
%             visc = alpha .*GR.dx^power.* abs(diff([bcVals, FF], 1, 2)./GR.dx).^2 ./ (0.5.*(abs(FF) + abs([bcVals, FF(:,1:end-1,:)]))) ;
            visc = alpha.* GR.dx.*abs(diff([bcVals, FF], 1, 2)).^2 ./ (0.5.*(abs(FF) + abs([bcVals, FF(:,1:end-1,:)]))) ;

        case 'E'
%             rhoSens = 1;%0.5.*(FF(:,:,indRho) + [FF(:,2:end, indRho), bcVals(:,:,indRho)]);
%             visc = alpha .*GR.dx^power.* abs(diff([FF, bcVals], 1, 2)./GR.dx).^2 ./ (0.5.*(abs(FF) + abs([FF(:,2:end, :), bcVals]))) ;
            visc = alpha .* GR.dx.*abs(diff([FF, bcVals], 1, 2)).^2 ./ (0.5.*(abs(FF) + abs([FF(:,2:end, :), bcVals]))) ;

        case 'N'
%             rhoSens = 1;%0.5.*(FF(:,:,indRho) + [FF(2:end,:,indRho); bcVals(:,:,indRho)]);
%             visc = alpha .*GR.dy^power.* abs(diff([FF; bcVals], 1, 1)./GR.dy).^2 ./ (0.5.*(abs(FF) + abs([FF(2:end,:,:); bcVals]))) ;
            visc = alpha.* GR.dy.*abs(diff([FF; bcVals], 1, 1)).^2 ./ (0.5.*(abs(FF) + abs([FF(2:end,:,:); bcVals]))) ;

        case 'S'
%             rhoSens = 1;%0.5.*(FF(:,:,indRho) + [bcVals(:,:,indRho); FF(1:end-1,:, indRho)]);
%             visc = alpha .*GR.dy^power.* abs(diff([bcVals; FF], 1, 1)./GR.dy).^2 ./ (0.5.*(abs(FF) + abs([bcVals; FF(1:end-1,:,:)]))) ;
            visc = alpha .*GR.dy.*abs(diff([bcVals; FF], 1, 1)).^2 ./ (0.5.*(abs(FF) + abs([bcVals; FF(1:end-1,:,:)]))) ;

    end
    testInd = isinf(visc) + isnan(visc);
    visc(logical(testInd)) = 0;
    visc = visc + epsMin;
end

end