%Stock market model with short-selling tax and endogenous shares: simulations 
%Last updated: Feb 14, 2024. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

clear, clc, %close all; 

%------------------
%Parameter values
%------------------
H = 100; 
r = 0.1; a = 1; 
betta = 2; 
dbar = 10; sigma = 1; Zbar = 0.1;  
pf = (dbar - a*sigma^2*Zbar)/r; %Fundamental price
Tax = 0.2; %Short-selling tax 
Tax_add = (1+r)*Tax/(a*sigma^2);
T = 100;  %no. of periods
No_Tax = 0; %Set No_Tax = 1 to simulate without short-selling tax (or set Tax = 0); 

%----------------
%Coding choices
%----------------
Fixed = 0; %Fixed  = 1: Pick fixed rather than time-varying (fitness-based) population shares. 
Naive = 0; %runs naive algorithm (Algo 1): starts from 1 non-buyer 

%-----------------
%Specify beliefs
%-----------------
%Disperse beliefs
b = zeros(H,1); C = b; g = b; g(ceil(H/2)+1:H) = 1 + linspace(0,0.4,H-ceil(H/2)); 
b(1:ceil(H/2)) = -0.2 + linspace(0,0.4,ceil(H/2)); 
C(1:ceil(H/2)) = 1-abs(b(1:ceil(H/2))); 

%No heterogeneity (two types)
%b = zeros(H,1); C = b; g = b; g(ceil(H/2)+1:H) = 1.2; C(1:ceil(H/2)) = 1-abs(b(1:ceil(H/2))); 

%Tax = 0.3, betta = 0.9,2,2.8,4
%Scenario 2 - betta = 2; dbar = 10; Tax = 0.30; p0 = pf + 1; H = 100;

%-------------------------------
%Initial values and predictors 
%-------------------------------
p0 = pf + 1; x0 = p0 - pf; xlag = p0 - pf; 
n_init = 1/H*ones(1,H);

%----------------------
%Preallocate matrices
%----------------------
U = NaN(H,1); Bind = NaN(T,1); Dlag2 = NaN(H,1); x = NaN(T,1); k_tot = x; k_tax = x; 
Beliefs = NaN(H,1);  AllNeg = zeros(T,1); AllZero = AllNeg; AllElse = AllNeg; 
Check1 = inf(T,1); Check11 = Check1;   %pop = NaN(H,T);

%--------------------------
%Generate dividend shocks 
%--------------------------
%Dividend shocks 
rng(1), sigma_d  = 0.01;
pd = makedist('Normal','mu',0,'sigma',sigma_d);  %Truncated normal distribution
pd_t = truncate(pd,-dbar,dbar); 
shock = random(pd_t,T,1);  
%shock = zeros(T,1); 

tic;  %Start stopclock

for t=1:T 
    
    %Beliefs = NaN(H,1);
    Divert = 0;
    
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

%--------------------------------------------       
%Obtain initial guess for no. non-buyers
%--------------------------------------------
   Demand_star = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*xstar)/(a*sigma^2);
   Demand_star(Demand_star<0 & Demand_star+Tax_add >=0) = 0;
   
   %Update k_init (uncomment to use)
   %if find(Demand_star==0) > 0
   %     zero = find(Demand_star==0);
   %     x_update = ( n_adj(1:min(zero)-1)*Beliefs_sort(1:min(zero)-1) + n_adj(max(zero)+1:end)*Beliefs_sort(max(zero)+1:end) + (1+r)*Tax*sum(n_adj(1:min(zero)-1))  - sum(n_adj(min(zero):max(zero)))*a*sigma^2*Zbar  ) / ( (1+r)*(1-sum(n_adj(min(zero):max(zero)))) );
   %else
   %     x_update = xstar + sum(n_adj(1:sum(Demand_star<0)))*Tax;
   %end
   %Demand_star = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*x_update)/(a*sigma^2);
   %Demand_star(Demand_star<0 & Demand_star+Tax_add >=0) = 0;
   
   k_init = sum(Demand_star<=0);

        if Naive == 1 
            k_init = 1;
        end

%--------------------------------------------------
%Measures of belief dispersion (Cases 2(i)-(ii))
%--------------------------------------------------

   disp0 = n_adj(k_init:end)*Beliefs_sort(k_init:end) - sum(n_adj(k_init:end))*Beliefs_sort(k_init);
   disp_hat = n_adj*Beliefs_sort - Beliefs_sort(k_init);
   sum_n = sum(n_adj(k_init:end));

   Stock_market_shorting_disp_insert

%-------------------------------------------------------------- 
%Find the equilibrium no. of short-sellers and zero positions
%--------------------------------------------------------------
  
        if AllZero(t)+AllNeg(t)==0 

            k_init(k_init==1) = 2;
            Stock_market_shorting_tax_cases_FAST
            %Stock_market_shorting_tax_cases

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

toc; %End stopclock


%Checks
max(Check1)
max(Check11)
sum(AllElse)
sum(Bind)

if sum(isnan(Check1)) + sum(isnan(Check11)) > 0
    disp('Market-clearing not satisfied')
end

if sum(Bind) == 0
    disp('Short-selling tax has no impact')
end

%---------------
%Plot figures
%--------------
%Time = 1:T;
%x_plot = [x0; x]; Time_plot = [0; Time']; 
%figure(1)
%hold on, subplot(2,2,2), plot(Time_plot,x_plot,'--k','LineWidth',1), hold on, 
%axis([-inf,inf,-inf,inf]), title('Scenario 2'), ylabel('Price deviation \it{x}'), xlabel('Time')


     


