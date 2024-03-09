%Stock market model with short-selling tax and endogenous shares: simulations 
%Insert for simulations using the file 'Stock_market_shorting_tax_BIFURC.m'
%Last updated: March 9, 2024. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

%------------------
%Parameter values
%------------------
%H = 100; 
%betta = 2;
r = 0.1; a = 1; 
dbar = 10; sigma = 1; Zbar = 0.1;  
pf = (dbar - a*sigma^2*Zbar)/r;  %Fundamental price
Tax = 0.1;  %Short-selling tax 
Tax_add = (1+r)*Tax/(a*sigma^2);

%----------------
%Coding choices
%----------------
%T = 100;  %no. of periods
No_Tax = 0; %Set No_Tax = 1 to simulate without short-selling tax (or set Tax = 0); 
Fixed = 0; %Fixed  = 1: Pick fixed rather than time-varying (fitness-based) population shares. 
Naive = 0; %runs naive algorithm (Algo 1): starts from 1 non-buyer 

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
        R = ones(length(Dlag2),1).*(x(t-1) + a*sigma^2*Zbar + shock(t-1) - (1+r)*x(t-2));

        if Bind(t-2) == 1
            if AllZero(t-2) == 1 
                Dlag2(Dlag2<0 & Dlag2+Tax_add >=0) = 0;
            elseif AllNeg(t) == 1 
                Dlag2(Dlag2+Tax_add<0) = Dlag2(Dlag2+Tax_add<0) + Tax_add;
            elseif AllElse(t) == 1 
            Dlag2(Dlag2<0 & Dlag2+Tax_add >=0) = 0;
            Dlag2(Dlag2+Tax_add<0) = Dlag2(Dlag2+Tax_add<0) + Tax_add;
            end
        R(Dlag2<0) = R(Dlag2<0) + (1+r)*Tax;
        end

        U = exp(betta*( R.*Dlag2 - C) );
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
    
%--------------------------------------------       
%Obtain initial guess for no. non-buyers
%--------------------------------------------
   Demand_star = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*xstar)/(a*sigma^2);
   Demand_star(Demand_star<0 & Demand_star+Tax_add >=0) = 0;

   %Stock_market_shorting_tax_k_update 
   %Uncomment to use
   
   k_sub = 0;    %Relevant when using 'Stock_market_shorting_tax_k_update'; else can set at 0
   k_init = max(sum(Demand_star<=0)-k_sub,1);

        if Naive == 1 
            k_init = 1;
        end

%--------------------------------------------------
%Measures of belief dispersion (Cases 2(i)-(ii))
%--------------------------------------------------

   sum_n = sum(n_adj(k_init:end));
   disp0 = n_adj(k_init:end)*Beliefs_sort(k_init:end) - sum_n*Beliefs_sort(k_init);
   disp_hat = n_adj*Beliefs_sort - Beliefs_sort(k_init);

   Stock_market_shorting_disp_insert

%-------------------------------------------------------------- 
%Find the equilibrium no. of short-sellers and zero positions
%--------------------------------------------------------------
  
        if AllZero(t)+AllNeg(t)==0 

            k_init(k_init==1) = 2;
            Stock_market_shorting_tax_cases_FAST
            %%Stock_market_shorting_tax_cases

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

    if sum(isnan(Demands))>0
       break   %break simulation if explosive path
    end

end




     


