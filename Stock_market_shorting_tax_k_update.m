%Stock_market_shorting_tax_k_update
%To be used with SIMS_FAST and SIMS_FAST_insert files
%Last updated: March 9, 2024. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

   if find(Demand_star==0) > 0
        zero = find(Demand_star==0);
        x_update = ( n_adj(1:min(zero)-1)*Beliefs_sort(1:min(zero)-1) + n_adj(max(zero)+1:end)*Beliefs_sort(max(zero)+1:end) + (1+r)*Tax*sum(n_adj(1:min(zero)-1))  - sum(n_adj(min(zero):max(zero)))*a*sigma^2*Zbar  ) / ( (1+r)*(1-sum(n_adj(min(zero):max(zero)))) );
   else
        x_update = xstar + sum(n_adj(1:sum(Demand_star<0)))*Tax;
   end
   Demand_star = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*x_update)/(a*sigma^2);
   Demand_star(Demand_star<0 & Demand_star+Tax_add >=0) = 0;