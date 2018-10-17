function Whittle_index_limited_PS = Whittle_indices_limited_PS(n,p,q, C0, C1, d, D)

[sum_Pi_n, sum_cost_n] = indices_limited_PS(n,p,q, C0, C1, d, D);
[sum_Pi_n_1, sum_cost_n_1] = indices_limited_PS(n-1,p,q, C0, C1, d, D);


sum_Pi_n-sum_Pi_n_1;
sum_cost_n - sum_cost_n_1;

Whittle_index_limited_PS = p*D -(sum_cost_n - sum_cost_n_1)/(sum_Pi_n - sum_Pi_n_1);

disp('Whittle index for limited PS')
disp(Whittle_index_limited_PS)


end
