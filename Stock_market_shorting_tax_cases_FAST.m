%Case of 2(iii) in Proosition 1: both short and zero position types 
%Last updated: Jan 2024. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

for k = k_init:length(Beliefs_sort)-1

    k_upper = k;

    for j=2:k_upper

         k_lower = j; 

         sum_n_tilde = sum(n_adj(k_lower:k_upper));  
         sum_n_tilde2 = sum(n_adj(1:k_lower-1));
         sum_Bel = n_adj(1:k_lower-1)*Beliefs_sort(1:k_lower-1) + n_adj(k_upper+1:end)*Beliefs_sort(k_upper+1:end);

         disp_tilde_init = sum_Bel - (1 - sum_n_tilde)*Beliefs_sort(k_upper);
         disp_tilde2_init = sum_Bel - (1 - sum_n_tilde)*Beliefs_sort(k_lower-1);

         disp_tilde = disp_tilde_init - (1-sum_n_tilde)*( Beliefs_sort(k_upper+1) - Beliefs_sort(k_upper) );
         disp_tilde2 = disp_tilde2_init - (1-sum_n_tilde)*( Beliefs_sort(k_lower) - Beliefs_sort(k_lower-1) );

         if max(disp_tilde,disp_tilde2-(1-sum_n_tilde)*(1+r)*Tax) <= a*sigma^2*Zbar - (1+r)*Tax*sum_n_tilde2  &&  min(disp_tilde_init, disp_tilde2_init - (1-sum_n_tilde)*(1+r)*Tax)  > a*sigma^2*Zbar - (1+r)*Tax*sum_n_tilde2
    
            k_tot(t) = k;  k_tax(t) = k_lower-1;
            x(t) = ( sum_Bel + (1+r)*Tax*sum_n_tilde2 - a*sigma^2*Zbar*sum_n_tilde ) / ( (1+r)*(1-sum_n_tilde) ); 
            AllElse(t) = 1;

         end

         if AllElse(t) == 1
             break
         end

    end

end



   
 