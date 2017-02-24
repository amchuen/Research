clear all;
close all;

% Test the thomas5 algorithm

% Create random pentadiagonal matrix
test_len = (1000);
max_val = 20;
test_mat = zeros(test_len);
loc = -1:1;
for i = 1:length(loc)
   test_mat = test_mat + diag(max_val * rand(1, test_len - abs(loc(i))), loc(i));
end

% Setup for thomas5 solver
ff = max_val * ones(test_len, 1);
% aa = [0; 0; diag(test_mat, -2)];
bb = [0; diag(test_mat, -1)];
cc = diag(test_mat);
dd = [diag(test_mat, 1); 0];
% ee = [diag(test_mat, 2); 0; 0];

% ff = [1, zeros(1, test_len-2), 1];
aa = zeros(1, test_len);
ee = zeros(1, test_len);
% bb = [0, -1*ones(1, test_len-1)];
% cc = 2*ones(1, test_len);
% dd = [-1*ones(1, test_len-1), 0];

disp(thomas5(aa, bb, cc, dd, ee, ff));

disp(thomas3(bb, cc, dd, ff));

% disp(test_mat\ff);