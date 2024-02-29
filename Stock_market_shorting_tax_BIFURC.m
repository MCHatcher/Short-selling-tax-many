%Stock market model with short-selling tax and endogenous shares: bifurcation diagrams 
%Last updated: Feb 15, 2024. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

clear, clc, %close all; 

%------------------
%Parameter values
%------------------
H = 20; 
window = 50;  %no. of terminal values to plot
T = 2350 + window;  %no. of periods
%n_betta = 100; 
n_betta = 82;
min_betta = 1; max_betta = 5;
n_sim = 20;

betta_vec = linspace(min_betta,max_betta,n_betta);

%Initial values and shocks
x_lower = -4;  x_upper = -0.01;
rng(5), x_init = x_lower + rand(n_sim,1)*(x_upper-x_lower);
shock = zeros(T,1);  %no dividend shocks

%Predictors
b = zeros(H,1); C = b; g = b; g(ceil(H/2)+1:H) = 1 + linspace(0,0.4,H-ceil(H/2)); 
b(1:ceil(H/2)) = -0.2 + linspace(0,0.4,ceil(H/2)); C(1:ceil(H/2)) = 1-abs(b(1:ceil(H/2))); 

%No heterogeneity (two types)
%b = zeros(H,1); C = b; g = b; g(ceil(H/2)+1:H) = 1.2; C(1:ceil(H/2)) = 1-abs(b(1:ceil(H/2))); 

%----------------------
%Preallocate matrices
%----------------------
brk = zeros(n_sim,1); C1 = NaN(n_sim,1); C11 = C1; C12 = C1; percent = NaN(n_betta,1); 
x_stack = NaN(window,n_sim); x_plot = NaN(window*n_sim,1); 

%------------------
%Run simulations
%------------------
%Coding choices
Fixed = 0;   %Fixed  = 1: Pick fixed rather than time-varying (fitness-based) population shares. 
Naive = 0;   %runs naive algorithm (Algo 1): starts from 1 non-buyer 

    for v = 1:length(betta_vec)

        betta = betta_vec(v);

            for s = 1:n_sim

                x0 = x_init(s); xlag = x0;

                n_init = 1/H*ones(1,H);

                Stock_market_shorting_tax_SIMS_FAST_insert
                %comment out betta, H and T

                %Store values for bifurc diagram
                x_stack(:,s) = x(end+1-window:end);  %Last 'window' observations

                if sum( isnan(x_stack(:,s)) ) == 0 && abs( range(x_stack(:,s)) ) < 5e-3
                    x_stack(:,s) = [NaN(window-1,1); x(end)];
                end

                %Check for non existence of an attractor
                r1 = 1 - isreal(x(end)); r2 = isnan(x(end));  r3 = isinf(x(end));

                if r1 + r2 + r3 > 0
                    brk(s) = 1;
                end

                %Accuracy measures
                C1(s) = max(Check1);
                C11(s) = max(Check11);
                C12(s) = max(Bind);

            end

        %Prepare to plot
        x_plot(:,v) = reshape(x_stack,1,[]);

        %Sims with no attractor
        percent(v) = 100*sum(brk)/n_sim;

        %Accuracy checks
        Check1_max = max(C1);
        Check11_max = max(C1);
        Bind_max = max(C12);

    end

%Final accuracy checks
max(Check1_max)
max(Check11_max)
max(Bind_max)

%Bifurcation plotter
figure(1)
subplot(1,2,1), 
hold on, xlabel('Intensity of choice \beta'), ylabel('Price deviation \it{x}')
axis([-inf,inf,-2.5,2.5])
plot(betta_vec,x_plot,'o','MarkerSize',2,'Color','k')  %'k' [0.5 0.5 0.5]

%for v=1:n_betta
%    plot(betta_vec(v),x_plot(:,v),'o','MarkerSize',2.2,'Color','k')
%end







     


