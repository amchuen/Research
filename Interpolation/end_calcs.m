function coeffs = end_calcs(xin, yin)
xleft = zeros(6,6);
yleft = zeros(6,1);
for i = 1:6
    yleft(i) = yin(i);
        for j = 1:6
            xleft(i,j) = xin(i)^(j-1);
        end
end

coeffs = xleft\yleft;
end