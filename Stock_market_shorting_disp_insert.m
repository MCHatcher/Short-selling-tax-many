%Stock_market_shorting_disp_insert 
%Insert for 'Stock_market_shorting_tax_SIMS_FAST.m' and 'Stock_market_shorting_tax_SIMS_FAST_insert.m'
%Last updated: March 9, 2024. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

       for m = k_init:length(Beliefs_sort)-1

           sum_n_init = sum_n;
           sum_n = sum_n - n_adj(m);
           Beliefs_diff = Beliefs_sort(m+1)- Beliefs_sort(m);
           Beliefs_diff2 = Beliefs_sort(1)- Beliefs_sort(m+1);

           disp_init = disp0;
           disp_hat_init = disp_hat;

           disp0 = disp0 - sum_n*Beliefs_diff;
           disp2 = disp0 - sum_n*Beliefs_diff2;
           %disp2 = n_adj(m+1:end)*Beliefs_sort(m+1:end) - sum(n_adj(m+1:end))*Beliefs_sort(1);
           disp_hat = disp_hat - Beliefs_diff;

            if max(disp0,disp2-(1+r)*Tax*sum_n) <= a*sigma^2*Zbar && disp_init > a*sigma^2*Zbar

                k_tot(t) = m;  k_tax(t) = 0;
                x(t) = ( n_adj(m+1:end)*Beliefs_sort(m+1:end) - ( 1-sum_n )*a*sigma^2*Zbar ) / ( (1+r)*sum_n );
                AllZero(t) = 1;
                break

            elseif disp_hat <=  a*sigma^2*Zbar - (1+r)*Tax*(1-sum_n)  &&  disp_hat_init - (1+r)*Tax > a*sigma^2*Zbar - (1+r)*Tax*(1-sum_n)   

                k_tot(t) = m;  k_tax(t) = m;
                x(t) = n_adj*Beliefs_sort/(1+r) + (1-sum_n)*Tax;
                AllNeg(t) = 1;
                break

            end

       end