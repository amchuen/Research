function [xout, yout, avec, bvec, cvec, dvec] = cubic_interp(xin, yin, res, k1, kn)
% This function will generate a cubic interplotaiton for a given set of
% data points. To accomplish this, the function will first generate a 

%% Numerical Interpolation?

% if strcmpi(type, 'num')
    % Use least squares approximation to determine the left side values
%     k_init = syseq_approx(xin(1:4), yin(1:4));
%     k_fin = syseq_approx(xin(end-3:end), yin(end-3:end));
[xout, yout, avec, bvec, cvec, dvec] = interp_fl1d(xin, yin, k1, kn, res);
% end

end

function kval = syseq_approx(xvals, yvals)

% use systems of equations to approximate first derivative for cubic
% functions... needs first four values

xleft = zeros(4,4);
for i = 1:4
   xleft(i,:) = [1, xvals(i), xvals(i)^2, xvals(i)^3];
end

coeff = xleft\yvals';
kval = coeff(2) + 2*coeff(3)*xvals(1) + 3*coeff(4)*xvals(1)^2;
end

function [xout, yout, avec, bvec, cvec, dvec] = interp_fl1d(xvals, yvals, k1, kn, res)
    % cubic interpolation with 1st derivative of FIRST and LAST POINT KNOWN
    
    % get each dx values
    hvals = diff(xvals);
    
    % get k values (first derivatives)
    alpha = [0, 1./hvals(1:end-1), 0];
    beta = [1, 2.0*(1./hvals(1:end-1) + 1./hvals(2:end)), 1];
    gamma = [0, 1./hvals(2:end), 0];
    delta = [k1, 3.*(((yvals(3:end) - yvals(2:end-1))./(hvals(2:end).^2)) + ((yvals(2:end-1) - yvals(1:end-2))./(hvals(1:end-1).^2))), kn];
    
    % use thomas algorithm to get the k vector using numerical equations
    kvals = thomas3(alpha, beta, gamma, delta);
    
    % initialize outputs
    xout = [];
    yout = [];
    avec = [];
    bvec = [];
    cvec = [];
    dvec = [];
    
    for i = 1:(length(xvals)-1)
        x_interp = linspace(xvals(i), xvals(i+1), res);

        aa3 = yvals(i);
        bb3 = kvals(i);
        cc3 = 3*(yvals(i+1) - yvals(i))/(hvals(i)^2) - (kvals(i+1)/hvals(i)) - 2*kvals(i)/hvals(i);
        dd3 = -2*(yvals(i+1) - yvals(i))/(hvals(i)^3) + (kvals(i+1) + kvals(i))/(hvals(i)^2);

        y_interp = aa3 + bb3.*(x_interp - xvals(i)) + cc3.*(x_interp-xvals(i)).^2 + dd3.*(x_interp-xvals(i)).^3;

        xout = [xout, x_interp];
        yout = [yout, y_interp];
        avec = [avec, aa3];
        bvec = [bvec, bb3];
        cvec = [cvec, cc3];
        dvec = [dvec, dd3];
    end
end