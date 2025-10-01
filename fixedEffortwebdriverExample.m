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
cd('DATA')
load('SimCons.mat')
cd('..')

rng(12390);
%for rep=1:1
rep=9999;
%% SIMULATIONS -------------------------------------------------------------------------
dt=1;
tspan=0:dt:150;

    %% SETUP
    k=randi(length(SimCons)); 
    web=SimCons(k).topo.web; %randomly choose a conserved web from the list used in this study
    spe=length(web);
    fish=SimCons(k).topo.fish;
    ext=SimCons(k).free.ext;
    B0=SimCons(k).free.B;
    
    T=SimCons(k).initial.Troph;
    [r,K,y,e,Bs,c,h,ax_ar,Z,po,Bext]=setup_default(web,fish);
    x=ax_ar.*(Z.^(-po.*(T-1)));
    
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
    tmp=false(nnz(fish'==1 & ext==0),1); 
    ind=randperm(nnz(fish'==1 & ext==0),nharvest);
    tmp(ind)=true;%use that index to mark that spot in the 0s vector
    harv(fish'==1 & ext==0)=tmp;%assign that to the harv list
    
    
    %% EFFORT LEVELS
    Effort=[0,1,2,3,4,5,6,7,8,9,10,12,15,20]; % list of fixed effort levels
    
    %% RUN
    E0=datasample(Effort,1); %randomly choose a Fixed Effort level, change this to look at specific levels
    sprintf('Simulation %d/%d Effort %d',k,length(SimCons),E0)
    X0=[B0,E0];

    %% old method
    %[t,X] = ode45(@(t,X) differential(t,X,x,y,r,K,e,h,c,Bs,web,harv,mu,co,ca,a,b,price,Bext),tspan,X0,options);

    n_steps = length(tspan);
    X = zeros(n_steps, length(X0));   % preallocate output
    X(1,:) = X0;                       % store initial state
    t = tspan(:);                      % store time points
    noise_scale1=0.05;
    noise_scale2=0.05;
    options=odeset('RelTol',10^-8,'AbsTol',10^-8);
    seasonality1=0.5;
    seasonality2=0.5;
    %disp(r);
    %disp(B);
    %disp(K);
    
    %disp(r.*(En-sum(B(basal))./K).*B);
    

for i = 2:n_steps
    noise=randn(size(B0));
    En1=noise'.*noise_scale1;
    En2=noise'.*noise_scale2;
    if mod(i,2)==0 %i is even
         Season1=1+seasonality1;
         Season2=1+seasonality2;
    else %i is odd
         Season1=1-seasonality1; 
         Season2=1-seasonality2; 
    end

    % code to extract the F values from the differential function each
    % timestep
    % Pack all parameters into a struct
    params = struct('x', x, 'y', y, 'r', r, 'K', K, 'e', e, 'h', h, 'c', c, ...
    'Bs', Bs, 'nicheweb', web, 'harvest', harv, 'mu', mu, ...
    'co', co, 'ca', ca, 'a', a, 'b', b, 'price', price, ...
    'Bext', Bext, 'En1', En1, 'En2', En2, 'Season1', Season1, 'Season2', Season2);

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

    B=X(:,1:spe);
    E=X(:,spe+1:end);
    B(B<Bext)=0;
    E(E<0)=0;
    X=[B,E];
    ext = X(end, :) == 0;
    count = sum(ext == 0);
    disp(count);
        

nSteps = length(t);
nSpecies = size(params.nicheweb, 1);
nStateVars = size(X, 2);

% Preallocate
F_series = zeros(nSpecies, nSpecies, nSteps);
J_series = zeros(nStateVars, nStateVars, nSteps);

for i = 1:nSteps
    xi = X(i,:)';      % state vector at time t(i), as column
    ti = t(i);         % current time

    % Extract biomass (you may need to adjust if B is not first in X)
    B = xi(1:nSpecies);  

    % Compute F and J
    F_series(:,:,i) = computeF(B, params);
    J_series(:,:,i) = computeJacobian(@differential, xi, ti, params);
end
%% PLOT RESULTS
%figure
%set(gcf,'color','w');
%Time series of all trophic species
%subplot(2,2,1)
%plot(X(:,1:30));
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
F_mean = squeeze(mean(F_series, 3));      % size: N × N
F_meanAbs = squeeze(mean(abs(F_series), 3));      % size: N × N
F_var  = squeeze(var(F_series, 0, 3));    % size: N × N
J_mean = squeeze(mean(J_series, 3));      % size: N × N
J_meanAbs = squeeze(mean(abs(J_series), 3));      % size: N × N
J_var  = squeeze(var(J_series, 0, 3));    % size: N × N



cd('Trials')
writematrix(T,"trophic_level"+rep+".csv");
writematrix(X,"foodweb_TS_FJ"+rep+".csv");
writematrix(X(:,fish),"foodweb_TS_fish"+rep+".csv");
writematrix(web,"foodweb_interactions_FJ"+rep+".csv");
writematrix(F_mean, "F_mean"+rep+".csv");
writematrix(F_var,  "F_var"+rep+".csv");
writematrix(J_mean, "J_mean"+rep+".csv");
writematrix(J_var,  "J_var"+rep+".csv");
writematrix(F_meanAbs,"F_meanAbs"+rep+".csv")
writematrix(J_meanAbs,"J_meanAbs"+rep+".csv")
tIndex=30;
writematrix(squeeze(F_series(:,:,tIndex)),"F_Inst_"+tIndex+"_9999.csv");
writematrix(squeeze(J_series(:,:,tIndex)),"J_Int_"+tIndex+"_9999.csv");
tIndex=40;
writematrix(squeeze(F_series(:,:,tIndex)),"F_Inst_"+tIndex+"_9999.csv");
writematrix(squeeze(J_series(:,:,tIndex)),"J_Int_"+tIndex+"_9999.csv");
tIndex=50;
writematrix(squeeze(F_series(:,:,tIndex)),"F_Inst_"+tIndex+"_9999.csv");
writematrix(squeeze(J_series(:,:,tIndex)),"J_Int_"+tIndex+"_9999.csv");
cd('..')

%end