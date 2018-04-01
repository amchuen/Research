function func = splinePts(coeffs, yVec)

func = zeros(size(yVec));

for i = 1:size(coeffs,2)
    func(i) = polyval(coeffs(:,i), yVec(i));
end


end