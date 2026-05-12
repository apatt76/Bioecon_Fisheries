%% constantwebdriverExample.m

%% Author --------------------------------------------------------------
% name: Valentin Cocco
% mail: valentin.cocco@ens.fr
% creation: 5-16-2018
% editor: Paul Glaum

%% Description ---------------------------------------------------------
% 1. Run the second stage of simulation with Fixed fishing effort on conserved foodwebs.

% Calls:
%   - webproperties.m
%   - setup_default.m
%   - differential.m

%% Updates -----------------------------------------------------------------
% when: 8-14-2019
% what: adaptation to updated analysis

%% LOAD DATA -------------------------------------------------------------------------------------
%This loads a file that contains the conserved webs from the fishing free films. 
cd('NewWebs')
load('Webs50.mat')
cd('..')

seed1=12398;
tmax=2000;
noise_scale1=0.05;
noise_scale2=0.05;

%rng(5);
for rep=1:50
new_seed=seed1;    
rng(new_seed);

    
%rep=100000;
%% SIMULATIONS -------------------------------------------------------------------------
dt=1;
tspan=0:dt:tmax;

    %% SETUP
    %k=randi(length(SimCons)); 
    k=rep;
    web=webs{k,1}; % choose a conserved web from the list used in this study
    spe=length(web);
    fish=webs{k,2};
    %ext=webs(k).free.ext;
    ext=zeros(1,30); % nobody is extinct yet
    B0=webs{k,3};
    
    T=webs{k,4};
    [r,K,y,e,Bs,c,h,ax_ar,Z,po,Bext]=setup_default(web,fish);
    % try averaging the fish and invertebrates:
    %ax_ar
    ax_ar(ax_ar~=0)=0.597;
    %y seems already same?
    


    x=ax_ar.*(Z.^(-po.*(T-1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %try to increase K
    K=K*2;
    % try to lower B_naught (denoted Bs here for B saturation)
    Bs=Bs*0.3;
    


    %mu=0, this keeps the Effort static at E0, the initial Effort value
    mu=0;
    
    co=1;
    ca=0.01;
    a=1;
    %b=0, this keeps the price static despite changing yield
    b=0;
    price='linear';
    nharvest=1;
    %randomly choose fish to harvest
    harv=false(spe,1);
    tmp=false(nnz( ext==0),1); 
    %ind=randperm(nnz(fish'==1 & ext==0),nharvest);
    ind=randperm(nnz(ext==0),nharvest); % harvest anything non extinct for now
    tmp(ind)=true;%use that index to mark that spot in the 0s vector
    harv( ext==0)=tmp;%assign that to the harv list
    
    
    %% EFFORT LEVELS
    %Effort=[0,1,2,3,4,5,6,7,8,9,10,12,15,20]; % list of fixed effort levels
    
    %% RUN
    %E0=datasample(Effort,1); %randomly choose a Fixed Effort level, change this to look at specific levels
    % set Effort to a constant
    E0 = 0; % Set a constant effort level for the simulation
    sprintf('Simulation %d/%d Effort %d',k,length(webs),E0)
    X0=[B0',E0];

    %% old method
    %[t,X] = ode45(@(t,X) differential(t,X,x,y,r,K,e,h,c,Bs,web,harv,mu,co,ca,a,b,price,Bext),tspan,X0,options);

    n_steps = length(tspan);
    X = zeros(n_steps, length(X0));   % preallocate output
    X(1,:) = X0;                       % store initial state
    t = tspan(:);                      % store time points

    options=odeset('RelTol',10^-8,'AbsTol',10^-8);
    seasonality1=0;
    seasonality2=0;
    %disp(r);
    %disp(B);
    %disp(K);
    
    %disp(r.*(En-sum(B(basal))./K).*B);
    


% while loop code
max_attempts = 10;  % max redraws
attempt = 0;
success = false;

while ~success && attempt < max_attempts
    attempt = attempt + 1;
    new_seed=seed1+attempt;
    rng(new_seed);
    % ---- 3a: create new random strengths ----
    strength1 = web;              % start with topology
    idx = (strength1 == 1);                       % existing links
    strength1=double(strength1);
    p1=0.1;
    p2=0.8;
    p3=0.1;
    rand_num1 = rand(sum(idx(:)),1);
    strength1(idx) = 0.5*(rand_num1 < p1) + ...
                 1.0*(rand_num1 >= p1 & rand_num1< p1+p2) + ...
                 1.5*(rand_num1 >= p1+p2);
    
    
    % ---- 3b: update params ---- 
    % include in 3c
 

    % set noise
    %eta=generate_OU_noise(30,n_steps,0,1,1,0);
    eta=randn(30,n_steps);
    eta(abs(eta)<1)=0.5*eta(abs(eta)<1); % make smallish noise even smaller
    
    
    % ---- 3c: integrate for next horizon ----
for i = 2:n_steps
    noise=eta(:, i);
    noise=noise.';
    %noise=zeros(1,30);
    En1=noise'.*noise_scale1+0.0;
    En2=noise'.*noise_scale2+0.0;
    if mod(i,2)==0 %i is even
         Season1=1+seasonality1;
         Season2=1+seasonality2;
    else %i is odd
         Season1=1-seasonality1; 
         Season2=1-seasonality2; 
    end


    if t(i)>10000 % for now strengths will not change, just run for 4000 TS
        strength=double(strength1);
    else
        strength=double(web);
    end

    % code to extract the F values from the differential function each
    % timestep
    % Pack all parameters into a struct
    params = struct('x', x, 'y', y, 'r', r, 'K', K, 'e', e, 'h', h, 'c', c, ...
    'Bs', Bs, 'nicheweb', web, 'harvest', harv, 'mu', mu, ...
    'co', co, 'ca', ca, 'a', a, 'b', b, 'price', price, ...
    'Bext', Bext, 'En1', En1, 'En2', En2, 'Season1', Season1, 'Season2', Season2, 'strength', strength);

    options = odeset('OutputFcn', @(t, y, flag) outputFvalues(t, y, flag, params));
    % Integrate from t(i-1) to t(i)
    [~, Xtemp] = ode45(@(t,X) differential(t,X,params),[t(i-1), t(i)], X(i-1,:)', options);
    
    % Take the final state and add noise (e.g., Gaussian noise)
    %disp(Xtemp(end,:)');
    
    
    %X_noisy = Xtemp(end,:)' + noise' .* noise_scale;  % tune noise_scale as needed
    %disp(X_noisy);
    % Save it
    X(i,:) = Xtemp(end,:)';
end
    
    % ---- 3d: check extinctions ----
    

    alive1 = X(1000,1:spe)>Bext;  % species alive at timepoint 1
    alive2 = X(tmax,1:spe)>Bext;  % species alive at timepoint 2

    % new extinctions = species that were alive at start but dead at end
    disp('here is the number of species alive');
    disp(sum(alive1));
    new_extinctions = alive1 & ~alive2;
    disp('here is the number of new extinctions');
    disp(sum(new_extinctions));
    % check if there were any new extinctions
    any_new_extinct = any(new_extinctions);
    
    %if ~any_new_extinct&&sum(alive1)>10
        success = true;
   %     fprintf('Success on attempt %d\n', attempt);
   %     disp('seed is');
   %     disp(new_seed);
   % else
   %     fprintf('Attempt did not pass muster %d, increasing seed\n', attempt);
   %     disp('seed is');
   %     disp(new_seed);
   % end
end



    B=X(:,1:spe);
    E=X(:,spe+1:end);
    B(B<Bext)=0;
    E(E<0)=0;
    X=[B,E];
    ext = X(end, :) == 0;
    count = sum(ext == 0);
    disp('number of species alive is');
    disp(count);
        

nSteps = length(t);
nSpecies = size(params.nicheweb, 1);
nStateVars = size(X, 2);

if success
    % Preallocate
%F_series = zeros(nSpecies, nSpecies, nSteps);
%J_series = zeros(nStateVars, nStateVars, nSteps);

%for i = (nSteps-2005):nSteps
%    xi = X(i,:)';      % state vector at time t(i), as column
%    ti = t(i);         % current time

    % Extract biomass (you may need to adjust if B is not first in X)
 %   B = xi(1:nSpecies);  

    % Compute F and J
 %   F_series(:,:,i) = computeF(B, params);
 %   J_series(:,:,i) = computeJacobian(@differential, xi, ti, params);
%end

%% PLOT RESULTS
figure
set(gcf,'color','w');
%Time series of all trophic species
%subplot(2,2,1)
plot(X(:,1:30));

%Time series of Effort (unchanging in the Fixed Effort simulation)
%subplot(2,2,2)
%plot(X(:,31));
%Time series of harvested species
%subplot(2,2,3)
%plot(X(:,harv))
%Time series of all labeled fish species
%subplot(2,2,4)
%plot(X(:,fish))



%disp(F_series);

%F_mean1 = squeeze(mean(F_series(:,:,2000:2200), 3));      % size: N × N
%F_meanAbs1 = squeeze(mean(abs(F_series(:,:,2000:2200)), 3));      % size: N × N
%F_var1  = squeeze(var(F_series, 0, 3));    % size: N × N
%J_mean1 = squeeze(mean(J_series(:,:,2000:2200), 3));      % size: N × N
%J_meanAbs1 = squeeze(mean(abs(J_series(:,:,2000:2200)), 3));      % size: N × N
%J_var1  = squeeze(var(J_series, 0, 3));    % size: N × N

%F_mean2 = squeeze(mean(F_series(:,:,3800:tmax), 3));      % size: N × N
%F_meanAbs2 = squeeze(mean(abs(F_series(:,:,3800:tmax)), 3));      % size: N × N
%F_var2  = squeeze(var(F_series, 0, 3));    % size: N × N
%J_mean2 = squeeze(mean(J_series(:,:,3800:tmax), 3));      % size: N × N
%J_meanAbs2 = squeeze(mean(abs(J_series(:,:,3800:tmax)), 3));      % size: N × N
%J_var2  = squeeze(var(J_series, 0, 3));    % size: N × N



cd('NewWebs')

exportgraphics(gca, "myPlot_connect50_"+rep+".png");

writematrix(X,"foodweb_TS_connect50_"+rep+"_"+new_seed+".csv");
%writematrix(T,"trophic_level"+rep+"_"+new_seed+".csv");
%writematrix(X(:,fish),"foodweb_TS_fish"+rep+"_"+new_seed+".csv");
%writematrix(web,"foodweb_interactions_FJ"+rep+"_"+new_seed+".csv");

%writematrix(strength1,"strength"+rep+"_"+new_seed+".csv")

%writematrix(J_mean1, "J_mean1_"+rep+"_"+new_seed+".csv");
%writematrix(J_var1,  "J_var1_"+rep+"_"+new_seed+".csv");
%writematrix(J_meanAbs1,"J_meanAbs1_"+rep+"_"+new_seed+".csv")
%writematrix(F_mean1, "F_mean1_"+rep+"_"+new_seed+".csv");
%writematrix(F_meanAbs1,"F_meanAbs1_"+rep+"_"+new_seed+".csv")
%{
writematrix(J_mean2, "J_mean2_"+rep+"_"+new_seed+".csv");
writematrix(J_var2,  "J_var2_"+rep+"_"+new_seed+".csv");
writematrix(J_meanAbs2,"J_meanAbs2_"+rep+"_"+new_seed+".csv");
writematrix(F_mean2, "F_mean2_"+rep+"_"+new_seed+".csv");
writematrix(F_meanAbs2,"F_meanAbs2_"+rep+"_"+new_seed+".csv")

writematrix(e,"e_"+rep+"_"+new_seed+".csv");
writematrix(ax_ar,"ax_ar"+rep+"_"+new_seed+".csv");
writematrix(y,"y_"+rep+"_"+new_seed+".csv");

[nrF, ncF, ntF] = size(F_series);
[nrJ, ncJ, ntJ] = size(J_series);

% output timeseries of J and F (optional)

F_out = zeros(ntF, nrF*ncF);
J_out = zeros(ntJ, nrJ*ncJ);


for tIndex = 1:ntF
    F_out(tIndex, :) = reshape(F_series(:,:,tIndex), 1, []);
end

for tIndex = 1:ntJ
    J_out(tIndex, :) = reshape(J_series(:,:,tIndex), 1, []);
end

writematrix(F_out, "F_Inst_" + rep +"_"+new_seed+ ".csv");
writematrix(J_out, "J_Inst_" + rep +"_"+new_seed+ ".csv");
%}
cd('..')
else
    % print message and go on to next one
    disp('web number');
    disp(rep);
end



end





