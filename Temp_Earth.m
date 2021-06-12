% Model the average temperature of the Earth.
% The model has N segments of athmosphere for which the temperature will be
% determined simultaneously with Earth's surface temperature.

clc
clear all
close all

%%

Atm.Height_atm = 50e3; % [m] Total height of athmosphere considered
Atm.Rair = 287.05 ; % [J/kg/K] Air specific gas constant = R/mu
Atm.N = 10; % Number of segments in the athmosphere
Atm.dh = Atm.Height_atm/Atm.N; % [m] Height of one segment of athmosphere

Atm.P0 = 101325; % [Pa] % Standard pressure at the sea level
Atm.I_Sun = 344; % [W/m^2] Flux coming from the Sun
Atm.alfa_sun = 0.33; % Reflected by the athmosphere/clouds etc
Atm.alfa_N = 0.04; % The earth’s surface (snow,ice,and so on) reflects a fraction alfa_N of what it receives
Atm.g = 9.81; % [m/s^2]


Atm.sigVis = 0.3*1e-4; % [m^-1] 
Atm.sigIR = 0.000188768; % [m^-1] % The larger the coeff the higher the absorbtion. 0 means no absorbtion

Atm.Bolt = 5.670374419184429453970e-8 ; %[W/m^2K^4]
Atm.Hydro = 113.6*0; %[W/m^2] Flux associated with hydrological scale. Moves heat from surface up to 10 km altitude (first 2 layers if total height=50 and N = 10)


T0 = 400*ones(1,Atm.N+1); % Initial values for iterations
lb = 10*ones(length(T0),1);
ub = (200+273)*ones(length(T0),1);


options = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point','UseParallel',true,'PlotFcn',{@optimplotx,@optimplotfval,@optimplotfirstorderopt},'StepTolerance',eps,'ConstraintTolerance',1e-14,'FiniteDifferenceType','central','MaxFunctionEvaluations',99999,'ScaleProblem',true,'OptimalityTolerance',1e-14); %'FiniteDifferenceStepSize',1e-15
[Topt, fval] = fmincon(@(T) error_eq_temp(T,Atm),T0,[],[],[],[],lb,ub,[],options);

% hybridopts = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point','UseParallel',true,'PlotFcn',{@optimplotx,@optimplotfval,@optimplotfirstorderopt},'StepTolerance',eps,'ConstraintTolerance',1e-14,'FiniteDifferenceType','forward','MaxFunctionEvaluations',99999,'ScaleProblem',true); %'FiniteDifferenceStepSize',1e-15
% options = optimoptions('particleswarm','SwarmSize', 50,'UseParallel',true,'HybridFcn',{@fmincon,hybridopts},'PlotFcn',{@pswplotbestf},'MaxStallIterations',300,'MinNeighborsFraction',0.35,'SocialAdjustmentWeight',1.2,'SelfAdjustmentWeight',1.2,'InertiaRange',[0.7,2.9]);
% [Topt, fval] = particleswarm(@(T) error_eq_temp(T,Atm),length(T0),lb,ub,options); 
% 





%%
for ii = 2:(Atm.N)
    T(ii) = {num2str(ii)};
end
T(1) = {'Upper Atm.'};
T(Atm.N+1) = {'Earth'};

X = categorical(T);
X = reordercats(X,T); 

figure
bar(X,Topt-273.15)
xlabel('Temperature in different layers')
ylabel('Temperature [\circC]')

%%


error_eq_temp(Topt,Atm)


%%
function ener_bala_err = error_eq_temp(T,Atm)

TN = T(end); % Temperature of the Earth
Trho = TN;
rho = @(h) Atm.P0/(Trho*Atm.Rair)*exp(-Atm.g*h/(Trho*Atm.Rair)); % Density function for the air given the temperature for the earth

% Average the density of air in a cell given its height
for ii = Atm.N:-1:1
    h0 = (Atm.N-ii)*Atm.dh;
    hf = (Atm.N-ii+1)*Atm.dh;
    domh = h0:1:hf;
    rho_ave(ii) = trapz(domh,rho(domh))/Atm.dh; % Height averaged Density of cell ii ;
end

% Visible light component going downwards from sun and component reflected from
% earth and going up

% Viz radiation
rhosum_above = cumsum(rho_ave);
rhosum_below = cumsum(rho_ave,'reverse');
F_viz_refl = Atm.alfa_N*(Atm.I_Sun*(1-Atm.alfa_sun)*exp(-Atm.sigVis*Atm.dh*rhosum_above(end)));

F_viz_transm = zeros(1,Atm.N+1);
F_viz_transm(1) =   Atm.I_Sun*(1-Atm.alfa_sun);
F_viz_transm_refl(Atm.N+1) =   F_viz_refl;

for ii = 2:(Atm.N+1) % Calculate transmitted Visible light and also the absorbed one. ii-1 is the cell, ii is the node before and after the cell
    F_viz_transm(ii) =   Atm.I_Sun*(1-Atm.alfa_sun)*exp(-Atm.sigVis*Atm.dh*rhosum_above(ii-1));
    F_viz_transm_refl(ii-1) =  F_viz_refl*exp(-Atm.sigVis*Atm.dh*rhosum_below(ii-1));
end
F_viz_absorbed_abov = -diff(F_viz_transm);
F_viz_absorbed_below =  diff(F_viz_transm_refl);


% IR radiation
Fabs_tot(Atm.N+1) = (1-Atm.alfa_N)*F_viz_transm(end); % The earth's absorbtion of visible radiation
for ii = 1:Atm.N %Loop through layers of air to calculate the total absorbed energy
    % What is absorbed (vis+IR) is emitted (IR)
    Fabs_tot(ii) = F_viz_absorbed_abov(ii)+F_viz_absorbed_below(ii); % The layers absorb visible light
    for jj = 1:Atm.N
        Fabs_tot(ii) = Fabs_tot(ii) + Fabs_ii_due_jj(T,ii,jj,Atm,rho_ave); % Add the absorbtion in IR due to other layers
    end
    Fabs_tot(ii) = Fabs_tot(ii) + TN^4*Atm.Bolt*( - exp(-Atm.sigIR*Atm.dh*sum(rho_ave(ii:end))) + exp(-Atm.sigIR*Atm.dh*sum(rho_ave((ii+1):end))));  % Add the contribution of earth. Earth emmits radiation only upwards
    
    Fabs_tot(Atm.N+1) = Fabs_tot(Atm.N+1) + 1/2*Atm.Bolt*T(ii)^4*exp(-Atm.sigIR*Atm.dh*sum(rho_ave((ii+1):end))); % Add to earth's absorbtion the contribution of each layer of air
end

F_tot_out_of_E = Atm.Bolt*T(end).^4*exp(-Atm.sigIR*Atm.dh*sum(rho_ave(:)));
for ii = 1:Atm.N
 F_tot_out_of_E = F_tot_out_of_E + 1/2*Atm.Bolt*T(ii).^4*exp(-Atm.sigIR*Atm.dh*sum(rho_ave(1:(ii-1))));
end
F_tot_out_of_E  = F_tot_out_of_E + F_viz_transm_refl(1); % Energy out of earth. Should be equal to 0.33*344 W/m^2 at internal thermal equilibrium

% Add hydro scale
Fabs_tot(end-1) = Fabs_tot(end-1) + Atm.Hydro/2; % Athmosphere close to the ground receives heat due to condensation
Fabs_tot(end-2) = Fabs_tot(end-2) + Atm.Hydro/2;
Fabs_tot(end)   = Fabs_tot(end)- Atm.Hydro; % Earth loses heat (through latent heat of vaporization)


ener_bala_err = Fabs_tot - Atm.Bolt*T.^4;
ener_bala_err = sum(ener_bala_err.^2);

end % error_eq_temp


function F_ir_abs = Fabs_ii_due_jj(T,ii,jj,Atm,rho_ave)
% Returns the flux absorbed in layer ii due to IR emission in layer jj.
% Earth is not considered an air layer
F_ir_E = 1/2*Atm.Bolt*T(jj)^4;
if ii<jj % Comes from below
    F_ir_abs =  F_ir_E*exp(-Atm.sigIR*Atm.dh*sum(rho_ave((ii+1):(jj-1))))  - F_ir_E*exp(-Atm.sigIR*Atm.dh*sum(rho_ave(ii:(jj-1))));
elseif ii>jj % Comes from above
    F_ir_abs =  F_ir_E*exp(-Atm.sigIR*Atm.dh*sum(rho_ave((jj+1):(ii-1))))  - F_ir_E*exp(-Atm.sigIR*Atm.dh*sum(rho_ave((jj+1):ii)));
else
    F_ir_abs = 0; % No absorbtion due to own emission
end
end % Fabs_ii_due_jj



