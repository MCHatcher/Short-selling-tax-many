%Stock market model with short-selling tax and endogenous shares: simulations 
%Last updated: Feb 14, 2024. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

%------------------
%Parameter values
%------------------
%H = 1000; 
%betta = 2; 
%T = 100;  %no. of periods
r = 0.1; a = 1; 
dbar = 10; sigma = 1; Zbar = 0.1;  
pf = (dbar - a*sigma^2*Zbar)/r; %Fundamental price
Tax = 0.1; %Short-selling tax 
Tax_add = (1+r)*Tax/(a*sigma^2);
No_Tax = 0; %Set No_Tax = 1 to simulate without short-selling tax (or set Tax = 0); 

%----------------------
%Preallocate matrices
%----------------------
U = NaN(H,1); Bind = NaN(T,1); Dlag2 = NaN(H,1); x = NaN(T,1); k_tot = x; k_tax = x; 
Beliefs = NaN(H,1);  AllNeg = zeros(T,1); AllZero = AllNeg; AllElse = AllNeg; 
Check1 = inf(T,1); Check11 = Check1;   %pop = NaN(H,T);

for t=1:T 
    
    %Beliefs = NaN(H,1); 
    
    if t==1
        Beliefs = b + g*x0;
        n = n_init;
    elseif t==2
        Beliefs = b + g*x(t-1);
        n = n_init; 
    elseif t>=3
        Beliefs = b + g*x(t-1);
        if t==3
            Dlag2 = (b + g*x0 + a*sigma^2*Zbar - (1+r)*x(t-2))/(a*sigma^2);
        else
            Dlag2 = (b + g*x(t-3) + a*sigma^2*Zbar - (1+r)*x(t-2))/(a*sigma^2);
        end
        if Bind(t-2) == 1
            if AllZero(t-2) == 1 
                Dlag2(Dlag2<0 & Dlag2+Tax_add >=0) = 0;
            elseif AllNeg(t) == 1 
                Dlag2(Dlag2+Tax_add<0) = Dlag2(Dlag2+Tax_add<0) + Tax_add;
            elseif AllElse(t) == 1 
            Dlag2(Dlag2<0 & Dlag2+Tax_add >=0) = 0;
            Dlag2(Dlag2+Tax_add<0) = Dlag2(Dlag2+Tax_add<0) + Tax_add;
            end
        end
        U = exp(betta*( (x(t-1) + a*sigma^2*Zbar + shock(t-1) - (1+r)*x(t-2))*Dlag2 - C) );
        n = transpose(U)/sum(U);

        %if Fixed==1
        %    n = n_init;
        %end

    end

    %pop(:,t) = n;

%------------------------------    
%Trial unconstrained solution
%------------------------------
xstar = n*Beliefs/(1+r);

if n*Beliefs - min(Beliefs) <= a*sigma^2*Zbar ||  No_Tax == 1    
    
        x(t) = xstar;   %Solution when taxes irrelevant or ignored 
        %Bind(t) = 0;

else 
        
        Bind(t) = 1;
        [Beliefs_sort,I] = sort(Beliefs);  
        n_adj = n(I);
        
%Sort beliefs when there are ties (uncomment to use, not essential)
        %if length(unique(Beliefs)) ~= H
        %    run Stock_market_shorting_sort_insert
        %end
    
        disp_stack = NaN(length(Beliefs_sort),1);
        disp_hat_stack = disp_stack; sum_stack = disp_stack;

%--------------------------------------------       
%Obtain initial guess for no. non-buyers
%--------------------------------------------
   Demand_star = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*xstar)/(a*sigma^2);
   k_init = sum(Demand_star<0); 

        if Naive==0
            
            if Iter == 1
                Stock_market_shorting_tax_iterations 
            else
                %Proceed
            end

        elseif Naive==1
            k_init = 1;
        end

    if k_init==0
        break
    end

%--------------------------------------------------
%Measures of belief dispersion (Cases 2(i)-(ii))
%--------------------------------------------------

       for m = k_init:length(Beliefs_sort)-1

           disp_init = n_adj(m:end)*Beliefs_sort(m:end) - sum(n_adj(m:end))*Beliefs_sort(m);
           disp = n_adj(m+1:end)*Beliefs_sort(m+1:end) - sum(n_adj(m+1:end))*Beliefs_sort(m+1);
           disp2 = n_adj(m+1:end)*Beliefs_sort(m+1:end) - sum(n_adj(m+1:end))*Beliefs_sort(1);

           disp_hat_init = n_adj*Beliefs_sort - Beliefs_sort(m);
           disp_hat = n_adj*Beliefs_sort - Beliefs_sort(m+1);
           sum_n = sum(n_adj(m+1:end));


            if max(disp,disp2-(1+r)*Tax*sum_n) <= a*sigma^2*Zbar && disp_init > a*sigma^2*Zbar

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

%-------------------------------------------------------------- 
%Find the equilibrium no. of short-sellers and zero positions
%--------------------------------------------------------------
  
    if AllZero(t)+AllNeg(t)==0 

        k_init(k_init==1) = 2;
        %Stock_market_shorting_tax_cases_FAST
        Stock_market_shorting_tax_cases

    end

end

%-----------------------
%Check market clearing
%-----------------------
Demands = (Beliefs + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);

if Bind(t) == 1

    Demands_adj = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    
    if AllZero(t) == 1 
        Demands(Demands<0 & Demands+Tax_add>=0)=0;
        Demands_adj(Demands_adj<0 & Demands_adj+Tax_add >=0)=0;
    elseif AllNeg(t) == 1
        Demands(Demands+Tax_add<0) = Demands(Demands+Tax_add<0) + Tax_add;
        Demands_adj(Demands_adj+Tax_add<0) = Demands_adj(Demands_adj+Tax_add<0) + Tax_add;
    elseif AllElse(t) == 1 
        Demands(Demands<0 & Demands+Tax_add>=0) = 0;
        Demands(Demands+Tax_add<0) = Demands(Demands+Tax_add<0) + Tax_add;
        Demands_adj(Demands_adj<0 & Demands_adj+Tax_add >=0) = 0;
        Demands_adj(Demands_adj+Tax_add<0) = Demands_adj(Demands_adj+Tax_add<0) + Tax_add;
    end

end

    Check1(t) = abs(n*Demands - Zbar);
    if Bind(t) == 1
        Check11(t) = abs(n_adj*Demands_adj - Zbar);
    else
        Check11(t) = Check1(t);
    end

end

if sum(isnan(Check1)) + sum(isnan(Check11)) > 0
    disp('Market-clearing not satisfied')
end

if sum(Bind) == 0
    disp('Short-selling tax has no impact')
end




     


