function [Uvals, CP] = flow_supersonic(flow, grid, yBu, yBl)

Uvals =  flow.U_inf * (1-cosd(flow.AOA)) .* ones(size(grid.XX));
vBtop = [0, diff(yBu)./grid.dx - flow.U_inf*sind(flow.AOA),  0];
vBbot = [0, diff(yBl)./grid.dx - flow.U_inf*sind(flow.AOA),  0];


%% Run Solver

beta2 = -(1 - flow.Ma^2);
vi = 2;

% Calculate U field using SLOR
for i = 3:(size( grid.XX, 2)) % loop through x direction
	% solve u_yy
    if (grid.xvals(i)>=0) % respect airfoil boundary condition
        AA = [0, (1/grid.dy^2) * ones(1, grid.ii_bot-1),... 
                0,...% top of airfoil
                (1/ grid.dy^2) * ones(1, length(grid.yvals)-grid.ii_top-1),0];
        BB = [1, (-beta2/(grid.dx^2) - 2/ grid.dy^2) * ones(1, grid.ii_bot-2),... 
                (-beta2/(grid.dx^2) - 1/ grid.dy^2),...% bottom of airfoil
                (-beta2/(grid.dx^2) - 1/ grid.dy^2),...% top of airfoil
                (-beta2/(grid.dx^2) - 2/ grid.dy^2) * ones(1, length(grid.yvals)-grid.ii_top-1),1];
        CC = [0, (1/ grid.dy^2) * ones(1, grid.ii_bot-2),...
                0,...% bottom of airfoil
                (1/ grid.dy^2) * ones(1, length(grid.yvals) - grid.ii_top),0];
        DD = [Uvals(1,i); beta2/( grid.dx^2).*(-2.*Uvals(2:( grid.ii_bot-1), i-1) + Uvals(2:( grid.ii_bot-1), i-2)); ...
                beta2/( grid.dx^2).*(-2.*Uvals(grid.ii_bot, i-1) + Uvals( grid.ii_bot, i-2)) - (vBbot(vi)-vBbot(vi-1))/( grid.dx* grid.dy);...
                beta2/( grid.dx^2).*(-2.*Uvals(grid.ii_top, i-1) + Uvals( grid.ii_top, i-2)) + (vBtop(vi)-vBtop(vi-1))/( grid.dx* grid.dy);...
                beta2/( grid.dx^2).*(-2.*Uvals((grid.ii_top+1):end-1, i-1) + Uvals(( grid.ii_top+1):end-1, i-2)); Uvals(end,i)];
        vi = vi + 1;

    else
        AA = [0, (1/grid.dy^2) * ones(1, length(grid.yvals)-2),0];
        BB = [1, (-beta2/(grid.dx^2) - 2/ grid.dy^2) * ones(1, length(grid.yvals)-2),1];
        CC = [0, (1/ grid.dy^2) * ones(1, length(grid.yvals)-2),0];
        DD = [Uvals(1,i); beta2/( grid.dx^2).*(-2.*Uvals(2:(end-1), i-1) + Uvals(2:(end-1), i-2)); Uvals(end,i)];
    end
%     AA = [0, (1/grid.dy^2).*ones(1, length(grid.yvals)-2), 0];
%     BB = [(beta/(grid.dx^2) - 1/(grid.dy^2)), (beta/(grid.dx^2) - 2/(grid.dy^2))*ones(1,length(grid.yvals)-2), 1];
%     CC = [(1/ grid.dy^2) * ones(1, length(grid.yvals) - 1),0];
%     DD = [(beta/(grid.dx^2)).*(2.*Uvals(1,i-1) - Uvals(1,i-2)) + (vBtop(i)-vBtop(i-1))/(grid.dx*grid.dy); (beta/(grid.dx^2)).*(2.*Uvals(2:(end-1),i-1) - Uvals(2:(end-1),i-2)); 1];

	Uvals(:,i) = thomas3(AA, BB, CC, DD);
    
%     CP = 1 - (Uvals.^2)./(flow.U_inf^2);
    CP = -2.*Uvals;
    
end

end
