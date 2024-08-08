function runMetabolicSwitching
clc
close all

% Choose reator type:
    reactorType = 'column'; % 'batch' or 'column'

% Choose model type:
    param.modelType = 'fba'; %'fba' or 'ann'

% Choose growth models:
    param.constraints = 'cybernetics'; % 'cybernetics' or 'kinetics'

% Choose plot options: 
    simulOptions.plot = 'on'; % 'on' or 'off;

% Determine the number of trials: 
 simulOptions.ntrial = 1; % any number >= 1

% basic parameter setting
param.kmax_lac = 22.1; 
param.kmax_pyr = 8.19; 
param.kmax_ace = 4.39; 
param.K_lac = 0.02; % mM
param.K_pyr = 0.02; % mM 
param.K_ace = 0.02; % mM
param.K_oxy = 0.006; % mM
param.k_d = 0.01; % cell death rate 


% parameters only for batch
param.kLa = 40; % 1/h
param.oxySat = 0.238; % mM

       
%% fba setting
switch param.modelType
    case 'ann'
        idx = [];

    case 'fba'
        % Input files
        model = readCbModel('iMR799_atp40.mat');
        mediaTbl = readtable('minimalMedia_noCarbon.xlsx');
        
        % Media condition
        idxExRxns = find(contains(model.rxns,'EX_'));
        model.lb(idxExRxns) = 0;
        idxRxnsInMedia = find(contains(model.rxns,mediaTbl.compounds));
        model.lb(idxRxnsInMedia) = -1000;
        
        % Rxn names and indices
        idx.rxnBio = find(strcmp(model.rxns,"bio1_biomass"));
        idx.rxnOxy = find(strcmp(model.rxns,"EX_cpd00007_e0"));
        idx.rxnAtp = find(strcmp(model.rxns,"ATPM_c0"));
        
        idx.rxnLac = find(strcmp(model.rxns,"EX_cpd00159_e0"));
        idx.rxnPyr = find(strcmp(model.rxns,"EX_cpd00020_e0"));
        idx.rxnAce = find(strcmp(model.rxns,"EX_cpd00029_e0"));
        
        idx.metAtp = find(contains(model.mets,"cpd00002_c0"));
        idx.metAdp = find(contains(model.mets,"cpd00008_c0"));
        idx.metH2o = find(contains(model.mets,"cpd00001_c0"));
        idx.metHplus = find(contains(model.mets,"cpd00067_c0"));

        % - parameters for FBA
        param.coeffAtp = 195.45;
        param.thresholdBioInLac = 0.6721; 
        param.thresholdPyrInLac = 0.6848;
        param.thresholdBioInPyr = 0.6837; 

        model.S(idx.metAtp,idx.rxnBio) = -param.coeffAtp;
        model.S(idx.metH2o,idx.rxnBio) = -param.coeffAtp;
        model.S(idx.metAdp,idx.rxnBio) = param.coeffAtp;
        model.S(idx.metHplus,idx.rxnBio) = param.coeffAtp;

        param.model = model;
        
end

switch reactorType
    case 'batch'
        T_elapsed = [];
        for isimul = 1:simulOptions.ntrial
            [t_elapsed nt] = goSimulation(reactorType,param,idx,simulOptions);
            T_elapsed = [T_elapsed;t_elapsed];
        end
        mean(T_elapsed)
        std(T_elapsed)
        nt
    
    case 'column'
%         for igrid = 20:10:100
        for igrid = 50
            simulOptions.ngrid = igrid;
            T_elapsed = [];
            for isimul = 1:simulOptions.ntrial
                [t_elapsed nt] = goSimulation(reactorType,param,idx,simulOptions);
                T_elapsed = [T_elapsed;t_elapsed];
            end
            igrid
            mean(T_elapsed)
            std(T_elapsed)
            nt
        end
end
%% 
function [t_elapsed nt] = goSimulation(reactorType,param,idx,simulOptions)
% parameters
odeOptions = odeset('NonNegative',1);

switch reactorType
    case 'batch'
        % lac, pyr, ace, oxy, biom
        y0 = [90 0 0 0.238 0.005];
        tspan = [0 35];

        tic
        [t,y] = ode45(@myOde_batch,tspan,y0,odeOptions,param,idx);
        % [t,y] = ode23s(@myOde,tspan,y0,odeOptions,param,idx);
        t_elapsed = toc;

        nt = length(t);

        if strcmp(simulOptions.plot,'on')
            figure
            subplot(3,1,1)
            plot(t,y(:,1:3),'linewidth',1.2)
            set(gca,'LineWidth',1.2,'Fontsize',12)
            ylim([0 100])
            xlim([tspan(1) tspan(end)])
            ylabel('Substrate [mM]')
            subplot(3,1,2)
            % semilogy(t,y(:,5),'linewidth',1.2)
            plot(t,y(:,5),'linewidth',1.2)
            set(gca,'LineWidth',1.2,'Fontsize',12)
            ylim([0 2.5])
            xlim([tspan(1) tspan(end)])
            ylabel('Biomass [g/L]')
            subplot(3,1,3)
            plot(t,y(:,4),'linewidth',1.2)
            set(gca,'LineWidth',1.2,'Fontsize',12)
            ylim([0 0.3])
            xlim([tspan(1) tspan(end)])
            ylabel('DO [mM]')
            xlabel('Time [hr]')
        end

    case 'column'
        param.columnHeight=0.3; % [m]
        param.nx = simulOptions.ngrid +1; % # of x_i = ngrid + 1
        param.delx = param.columnHeight/(param.nx-1); % delta x 
%         param.v=0.05; % velocity [m/d]
        param.v=0.01; % velocity [m/d]
        param.D=0.015; % diffusivity [m^2/d] 
        param.porosity=0.42; % porosity (Liu et al., 2017; Table S2)
        % para.porosity=1; % porosity
        
        % lac, pyr, ace, oxy, biom
        param.y_feed = [90 0 0 0.238 0];

        y0_lac = zeros(1,param.nx);
        y0_pyr = zeros(1,param.nx);
        y0_ace = zeros(1,param.nx);
        y0_oxy = param.oxySat*ones(1,param.nx);
        y0_bio = 0.005*ones(1,param.nx);

        y0 = [y0_lac y0_pyr y0_ace y0_oxy y0_bio];

        tspan = [0,5];
%         tspan = 0:1:10;
        
        tic
        [t,y] = ode45(@myOde_column,tspan,y0,odeOptions,param,idx);
        t_elapsed = toc;

        nt = length(t);

        if strcmp(simulOptions.plot,'on')
            nt = length(t);
            deltt = floor(nt/5);
            xx=0:param.delx:param.columnHeight;
            figure
            for it=1:deltt:nt
                for j=1:5
                    subplot(5,1,j)
                    plot(xx,y(it,(j-1)*param.nx+1:j*param.nx))
                    drawnow 
                    hold on
                end
            end
        end

end

%%-------------------------------------------------------------------------
function dydt = myOde_batch(t,y,param,idx)
ny = length(y);
dydt = zeros(ny,1);

% idxNeg = y<0;
% y(idxNeg)=0;

input_lac = zeros(2,1); % uptake rates of lac and oxy (minus sign)
input_pyr = zeros(2,1); % uptake rates of pyr and oxy (minus sign)
input_ace = zeros(2,1); % uptake rates of ace and oxy (minus sign)

% convert y to physical vars
lac = y(1);
pyr = y(2);
ace = y(3);
oxy = y(4);
bio = y(5);

modelType = param.modelType;
constraints = param.constraints;

kmax_lac = param.kmax_lac; 
kmax_pyr = param.kmax_pyr; 
kmax_ace = param.kmax_ace; 
K_lac = param.K_lac; 
K_pyr = param.K_pyr; 
K_ace = param.K_ace; 
K_oxy = param.K_oxy; 
k_d = param.k_d;
kLa = param.kLa;
oxySat = param.oxySat;

% inputs
input_lac(1) = kmax_lac*lac/(K_lac + lac); % lactate consumption 
input_lac(2) = kmax_lac*oxy/(K_oxy + oxy); % oxygen consumption (together with lactate) 
input_pyr(1) = kmax_pyr*pyr/(K_pyr + pyr); % pyruvate consumption 
input_pyr(2) = kmax_pyr*oxy/(K_oxy + oxy); % oxygen consumption (together with pyruvate)
input_ace(1) = kmax_ace*ace/(K_ace + ace); % acetate consumption 
input_ace(2) = kmax_ace*oxy/(K_oxy + oxy); % oxygen consumption (together with acetate)

% cybernetic variables 
switch constraints 
    case 'kinetics'
        u_lac = 1;
        u_pyr = 1;
        u_ace = 1;
    case 'cybernetics'
        % cybernetic variables (obj: carbon influx)
        den = 3*input_lac(1) + 3*input_pyr(1) + 2*input_ace(1);
        u_lac = 3*input_lac(1)/den;
        u_pyr = 3*input_pyr(1)/den;
        u_ace = 2*input_ace(1)/den;
end

% adjusted inputs
input_lac = input_lac*u_lac;
input_pyr = input_pyr*u_pyr;
input_ace = input_ace*u_ace;


switch modelType
    case 'fba'

        %----
        if input_lac(1)<1e-8
            input_lac(1)=0;
        end
        if input_pyr(1)<1e-8
            input_pyr(1)=0;
        end
        if input_ace(1)<1e-8
            input_ace(1)=0;
        end


        model = param.model;

        thresholdBioInLac = param.thresholdBioInLac; 
        thresholdPyrInLac = param.thresholdPyrInLac;
        thresholdBioInPyr = param.thresholdBioInPyr; 

        uptake.lac = input_lac(1);
        uptake.pyr = input_pyr(1);
        uptake.ace = input_ace(1);
        
        % growth on lactate 
        uptake.oxy = input_lac(2);
        x_idx_lac = goLac(model,idx,uptake,thresholdBioInLac,thresholdPyrInLac);
        
        % growth on pyruvate
        uptake.oxy = input_pyr(2);
        x_idx_pyr = goPyr(model,idx,uptake,thresholdBioInPyr);
        
        % growth on acetate 
        uptake.oxy = input_ace(2);
        x_idx_ace = goAce(model,idx,uptake);

        rlac_lac = - input_lac(1); % uptake
        roxy_lac = x_idx_lac(4); % minus sign (=uptake)
        rbio_lac = x_idx_lac(5);
        rpyr_lac = x_idx_lac(2);
        race_lac = x_idx_lac(3);

        rpyr_pyr = - input_pyr(1); % uptake
        roxy_pyr = x_idx_pyr(4); % minus sign (=uptake)
        rbio_pyr = x_idx_pyr(5);
        race_pyr = x_idx_pyr(3);
        
        race_ace = - input_ace(1); % uptake
        roxy_ace = x_idx_ace(4); % minus sign (=uptake)
        rbio_ace = x_idx_ace(5);
        
    case 'ann'
        output_lac = myAnn_lac_all(input_lac); % output: oxy(minus), bio, pyr, ace
        output_pyr = myAnn_pyr_all(input_pyr); % output: oxy(minus), bio, ace
        output_ace = myAnn_ace_all(input_ace); % output: oxy(minus), bio

        rlac_lac = - input_lac(1); % uptake
        roxy_lac = output_lac(1); % minus sign (=uptake)
        rbio_lac = output_lac(2);
        rpyr_lac = output_lac(3);
        race_lac = output_lac(4);

        rpyr_pyr = - input_pyr(1); % uptake
        roxy_pyr = output_pyr(1); % minus sign (=uptake)
        rbio_pyr = output_pyr(2);
        race_pyr = output_pyr(3);
        
        race_ace = - input_ace(1); % uptake
        roxy_ace = output_ace(1); % minus sign (=uptake)
        rbio_ace = output_ace(2);
        
end

% ODEs: lac, pyr, ace, oxy, bio
dydt(1) = rlac_lac*bio; % lac consumption 
dydt(2) = (rpyr_lac + rpyr_pyr)*bio; % pyr consumption or production 
dydt(3) = (race_lac + race_pyr + race_ace)*bio; % ace consumption or production 
dydt(4) = (roxy_lac + roxy_pyr + roxy_ace)*bio + kLa*(oxySat - oxy); % oxy consumption 
dydt(5) = (rbio_lac + rbio_pyr + rbio_ace)*bio - k_d*bio; % bio production 

%%-------------------------------------------------------------------------
function dydt = myOde_column(t,y,param,idx)
ny = length(y); dydt = zeros(ny,1);
% idxNeg = y<0; y(idxNeg)=0;

columnHeight = param.columnHeight;
nx = param.nx;
delx = param.delx;
v = param.v;
D = param.D;
porosity= param.porosity;
y_feed = param.y_feed; % lac, pyr, ace, oxy, biom

modelType = param.modelType;
constraints = param.constraints;

kmax_lac = param.kmax_lac; 
kmax_pyr = param.kmax_pyr; 
kmax_ace = param.kmax_ace; 
K_lac = param.K_lac; 
K_pyr = param.K_pyr; 
K_ace = param.K_ace; 
k_d = param.k_d;
K_oxy = param.K_oxy; 
kLa = param.kLa;
oxySat = param.oxySat;

input_lac = zeros(2,1); % uptake rates of lac and oxy (minus sign)
input_pyr = zeros(2,1); % uptake rates of pyr and oxy (minus sign)
input_ace = zeros(2,1); % uptake rates of ace and oxy (minus sign)

for ix = 1:nx

    % convert y to physical vars
    lac = y(ix);
    pyr = y(nx+ix);
    ace = y(2*nx+ix);
    oxy = y(3*nx+ix);
    bio = y(4*nx+ix);

    if ix==1 % Dackwerts B.C.
        lac_left = y(ix) - v/(porosity*D)*(y(ix)-y_feed(1));
        pyr_left = y(nx+ix) - v/(porosity*D)*(y(nx+ix)-y_feed(2));
        ace_left = y(2*nx+ix) - v/(porosity*D)*(y(2*nx+ix)-y_feed(3));
        oxy_left = y(3*nx+ix) - v/(porosity*D)*(y(3*nx+ix)-y_feed(4));
        bio_left = y(4*nx+ix) - v/(porosity*D)*(y(4*nx+ix)-y_feed(5));
    else
        lac_left = y(ix-1);
        pyr_left = y(nx+ix-1);
        ace_left = y(2*nx+ix-1);
        oxy_left = y(3*nx+ix-1);
        bio_left = y(4*nx+ix-1);
    end

    if ix==nx % Neumann (no flux) B.C. 
        lac_right = y(ix-1);
        pyr_right = y(nx+ix-1);
        ace_right = y(2*nx+ix-1);
        oxy_right = y(3*nx+ix-1);
        bio_right = y(4*nx+ix-1);        
    else
        lac_right = y(ix+1);
        pyr_right = y(nx+ix+1);
        ace_right = y(2*nx+ix+1);
        oxy_right = y(3*nx+ix+1);
        bio_right = y(4*nx+ix+1);
    end
    
    % inputs
    input_lac(1) = kmax_lac*lac/(K_lac + lac);
    input_lac(2) = kmax_lac*oxy/(K_oxy + oxy);
    input_pyr(1) = kmax_pyr*pyr/(K_pyr + pyr);
    input_pyr(2) = kmax_pyr*oxy/(K_oxy + oxy);
    input_ace(1) = kmax_ace*ace/(K_ace + ace);
    input_ace(2) = kmax_ace*oxy/(K_oxy + oxy);
    
    switch constraints 
        case 'kinetics'
            u_lac = 1;
            u_pyr = 1;
            u_ace = 1;
        case 'cybernetics'
            % cybernetic variables (obj: carbon influx)
            den = 3*input_lac(1) + 3*input_pyr(1) + 2*input_ace(1);
            if den < eps
                u_lac = 0;
                u_pyr = 0;
                u_ace = 0;
            else
                u_lac = 3*input_lac(1)/den;
                u_pyr = 3*input_pyr(1)/den;
                u_ace = 2*input_ace(1)/den;
            end
    end
    
    % adjusted inputs
    input_lac = input_lac*u_lac;
    input_pyr = input_pyr*u_pyr;
    input_ace = input_ace*u_ace;    
    
    switch modelType
        case 'fba'
    
            %----
            if input_lac(1)<1e-8
                input_lac(1)=0;
            end
            if input_pyr(1)<1e-8
                input_pyr(1)=0;
            end
            if input_ace(1)<1e-8
                input_ace(1)=0;
            end    
    
            model = param.model;
    
            thresholdBioInLac = param.thresholdBioInLac; 
            thresholdPyrInLac = param.thresholdPyrInLac;
            thresholdBioInPyr = param.thresholdBioInPyr; 
    
            uptake.lac = input_lac(1);
            uptake.pyr = input_pyr(1);
            uptake.ace = input_ace(1);
            
            % growth on lactate 
            uptake.oxy = input_lac(2);
            x_idx_lac = goLac(model,idx,uptake,thresholdBioInLac,thresholdPyrInLac);
            
            % growth on pyruvate
            uptake.oxy = input_pyr(2);
            x_idx_pyr = goPyr(model,idx,uptake,thresholdBioInPyr);
            
            % growth on acetate 
            uptake.oxy = input_ace(2);
            x_idx_ace = goAce(model,idx,uptake);
    
            rlac_lac = - input_lac(1); % uptake
            roxy_lac = x_idx_lac(4); % minus sign (=uptake)
            rbio_lac = x_idx_lac(5);
            rpyr_lac = x_idx_lac(2);
            race_lac = x_idx_lac(3);
    
            rpyr_pyr = - input_pyr(1); % uptake
            roxy_pyr = x_idx_pyr(4); % minus sign (=uptake)
            rbio_pyr = x_idx_pyr(5);
            race_pyr = x_idx_pyr(3);
            
            race_ace = - input_ace(1); % uptake
            roxy_ace = x_idx_ace(4); % minus sign (=uptake)
            rbio_ace = x_idx_ace(5);
            
        case 'ann'
            output_lac = myAnn_lac_all(input_lac); % output: oxy(minus), bio, pyr, ace
            output_pyr = myAnn_pyr_all(input_pyr); % output: oxy(minus), bio, ace
            output_ace = myAnn_ace_all(input_ace); % output: oxy(minus), bio
    
            rlac_lac = - input_lac(1); % uptake
            roxy_lac = output_lac(1); % minus sign (=uptake)
            rbio_lac = output_lac(2);
            rpyr_lac = output_lac(3);
            race_lac = output_lac(4);

            rpyr_pyr = - input_pyr(1); % uptake
            roxy_pyr = output_pyr(1); % minus sign (=uptake)
            rbio_pyr = output_pyr(2);
            race_pyr = output_pyr(3);
            
            race_ace = - input_ace(1); % uptake
            roxy_ace = output_ace(1); % minus sign (=uptake)
            rbio_ace = output_ace(2);            
    end
    
    % ODEs: lac, pyr, ace, oxy, bio
    diff_lac = D*(lac_right - 2*lac + lac_left)/delx^2;
    diff_pyr = D*(pyr_right - 2*pyr + pyr_left)/delx^2;
    diff_ace = D*(ace_right - 2*ace + ace_left)/delx^2;
    diff_oxy = D*(oxy_right - 2*oxy + oxy_left)/delx^2;
%     diff_bio = D*(bio_right - 2*bio + bio_left)/delx^2;
    diff_bio = 0;
    
    conv_lac = (v/porosity)*(lac - lac_left)/delx;
    conv_pyr = (v/porosity)*(pyr - pyr_left)/delx;
    conv_ace = (v/porosity)*(ace - ace_left)/delx;
    conv_oxy = (v/porosity)*(oxy - oxy_left)/delx;
%     conv_bio = (v/porosity)*(bio - bio_left)/delx;
    conv_bio = 0;

    sour_lac = rlac_lac*bio; % lac consumption
    sour_pyr = (rpyr_lac + rpyr_pyr)*bio; % pyr consumption or production
    sour_ace = (race_lac + race_pyr + race_ace)*bio; % ace consumption or production
%     sour_oxy = (roxy_lac + roxy_pyr + roxy_ace)*bio + kLa*(oxySat - oxy); % oxy consumption
    sour_oxy = (roxy_lac + roxy_pyr + roxy_ace)*bio; % oxy consumption without kLa term
    sour_bio = (rbio_lac + rbio_pyr + rbio_ace)*bio - k_d*bio; % bio production
    
    dydt(ix) = diff_lac - conv_lac + sour_lac;
    dydt(nx+ix) = diff_pyr - conv_pyr + sour_pyr;
    dydt(2*nx+ix) = diff_ace - conv_ace + sour_ace;
    dydt(3*nx+ix) = diff_oxy - conv_oxy + sour_oxy;
    dydt(4*nx+ix) = diff_bio - conv_bio + sour_bio;

end

%% 
function x_idx = goLac(model,idx,uptake,thresholdBio,thresholdPyr)

[nmets,nrxns] = size(model.S);

model.lb(idx.rxnLac) = -uptake.lac; % lactate
model.lb(idx.rxnOxy) = -uptake.oxy;

% maximize biomass 
f = zeros(nrxns,1);
f(idx.rxnBio) = -1; % biomass maximization

Aineq = []; bineq = [];
Aeq = model.S; beq = zeros(nmets,1);
[x,fval,exitflag] = cplexlp(f,Aineq,bineq,Aeq,beq,model.lb,model.ub);
if exitflag <= 0
%     exitflag
%     warning('incorrect solution...')
    bioMax = 0;
else
    bioMax = x(idx.rxnBio);
end

% maximize pyruvate
f = zeros(nrxns,1); f(idx.rxnPyr) = -1; % pyruvate maximization
model.lb(idx.rxnBio) = thresholdBio*bioMax; % constrain biomass <---

[x,fval,exitflag] = cplexlp(f,Aineq,bineq,Aeq,beq,model.lb,model.ub);
if exitflag <= 0
%     exitflag
%     warning('incorrect solution...')
    pyrMax = 0;
else
    pyrMax = x(idx.rxnPyr);
end

% maximize acetate
f = zeros(nrxns,1); f(idx.rxnAce) = -1; % pyruvate maximization
model.lb(idx.rxnPyr) = thresholdPyr*pyrMax; % constrain pyruvate production <---

[x,fval,exitflag] = cplexlp(f,Aineq,bineq,Aeq,beq,model.lb,model.ub);
if exitflag <= 0
%     exitflag
%     warning('incorrect solution...')
    x = zeros(length(x),1);
end
x_idx = [x(idx.rxnLac) x(idx.rxnPyr) x(idx.rxnAce) x(idx.rxnOxy) x(idx.rxnBio)];

%%
function x_idx = goPyr(model,idx,uptake,thresholdBio)

[nmets,nrxns] = size(model.S);

model.lb(idx.rxnPyr) = -uptake.pyr; % pyruvate
model.lb(idx.rxnOxy) = -uptake.oxy;

% maximize biomass 
f = zeros(nrxns,1);
f(idx.rxnBio) = -1; % biomass maximization

Aineq = []; bineq = [];
Aeq = model.S; beq = zeros(nmets,1);
[x,fval,exitflag] = cplexlp(f,Aineq,bineq,Aeq,beq,model.lb,model.ub);
if exitflag <= 0
%     exitflag
%     warning('incorrect solution...')
end
bioMax = x(idx.rxnBio);

% maximize acetate
f = zeros(nrxns,1); f(idx.rxnAce) = -1; % pyruvate maximization
model.lb(idx.rxnBio) = thresholdBio*bioMax; % constrain biomass <---

[x,fval,exitflag] = cplexlp(f,Aineq,bineq,Aeq,beq,model.lb,model.ub);
if exitflag <= 0
%     exitflag
%     warning('incorrect solution...')
    x = zeros(length(x),1);
end
x_idx = [x(idx.rxnLac) x(idx.rxnPyr) x(idx.rxnAce) x(idx.rxnOxy) x(idx.rxnBio)];

%%
function x_idx = goAce(model,idx,uptake)

[nmets,nrxns] = size(model.S);

model.lb(idx.rxnAce) = -uptake.ace; % pyruvate
model.lb(idx.rxnOxy) = -uptake.oxy;

% maximize biomass
f = zeros(nrxns,1);
f(idx.rxnBio) = -1; % biomass maximization

Aineq = []; bineq = [];
Aeq = model.S; beq = zeros(nmets,1);
[x,fval,exitflag] = cplexlp(f,Aineq,bineq,Aeq,beq,model.lb,model.ub);
if exitflag <= 0
%     exitflag
%     warning('incorrect solution...')
    x = zeros(length(x),1);
end
x_idx = [x(idx.rxnLac) x(idx.rxnPyr) x(idx.rxnAce) x(idx.rxnOxy) x(idx.rxnBio)];