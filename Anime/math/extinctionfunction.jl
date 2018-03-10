using Distributions

numpreds = collect(0:0.1:10);
lpreds = length(numpreds);
probext = Array{Float64}(lpreds);
P = 0.1;
r = 0.3;
lr = 10000;
for i=1:lpreds
    mu = 1 - ((numpreds[i]*P)/r);
    ndist = Normal(mu,1);
    draws = rand(ndist,lr);
    probext[i] = length(find(x->x<0,draws))/lr;
end

R"plot($numpreds,$probext,ylim=c(0,1))"

