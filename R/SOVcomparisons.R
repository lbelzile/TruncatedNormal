#precisionSOV vs covarianceSOV
#estimate numerically the probability that X≥0d for X∼N(0d,0.5Id+0.51d1⊤d), which is known and equal (d+1)^−1.

d = 10;
l = rep(0,d);
u = rep(Inf,d);
sigma = 0.5 * (diag(d) + matrix(1, d, d));
sigmaInv = solve(sigma);

n = 10000; #replications

est1 = 0;

for(i in (1:n)){
  est1 = est1 + precisionSOV(sigmaInv, l, u);
}
est1 = est1/n;

est2 = 0;

for(i in (1:n)){
  est2 = est2 + covarianceSOV(sigma, l, u);
}

est2 = est2/n;

print("precisionSOV error: ");
print(abs(est1 - 1/(d+1))*(d+1));

print("covarianceSOV error: ");
print(abs(est2 - 1/(d+1))*(d+1));
