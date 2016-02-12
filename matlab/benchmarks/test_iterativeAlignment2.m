% Simulate contour length increments
clear r
trueCidx = [];
lcSD=2;

for i=1:20
    r{i}.L = [ 30 + lcSD.*randn(1,1) 70 + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1)]+ -10 + (30--10)*rand(1,1);
    trueCidx(i) = 1;
end
cidx = ones(size(r));
perm = ones(size(r));
iterations = 1;
group = ones(size(r));
iterativeAlignment2(r,ones(length(r)),ones(length(r)),1,ones(length(r)))
