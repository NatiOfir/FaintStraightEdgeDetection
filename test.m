tic;
parfor i = 1:10^3
    A = rand(1000);
    A = inv(A);
end
toc;