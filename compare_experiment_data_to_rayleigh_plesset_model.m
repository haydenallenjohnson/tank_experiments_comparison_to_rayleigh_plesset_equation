%% load experiment data
load('experiment_data_20240604.mat');
experiment_data = struct;

experiment_data(2).scale_reading = data1(:,1)*1e-3; % m
experiment_data(2).rms_volts_reference = data1(:,2)*1e-3; % V
experiment_data(2).rms_volts_signal = data1(:,3)*1e-3; % V

experiment_data(1).scale_reading = data2(:,1)*1e-3; % m
experiment_data(1).rms_volts_reference = data2(:,2)*1e-3; % V
experiment_data(1).rms_volts_signal = data2(:,3)*1e-3; % V

% fit over spherical spreading regime to determine scale offset
xmin = 2e-3;
xmax = 10e-3;

for i = 1:length(experiment_data)
    experiment_data(i).signal_ratio = experiment_data(i).rms_volts_signal./experiment_data(i).rms_volts_reference;
    fit_ind = experiment_data(i).scale_reading >= xmin & experiment_data(i).scale_reading <= xmax;
    coef = polyfit(experiment_data(i).scale_reading(fit_ind), 1./experiment_data(i).signal_ratio(fit_ind),1);
    m = coef(1);
    b = coef(2);
    A = 1/m;
    experiment_data(i).scale_offset = A*b;
    experiment_data(i).distance = experiment_data(i).scale_reading + experiment_data(i).scale_offset;
    experiment_data(i).signal_ratio_times_distance = mean(experiment_data(i).distance(fit_ind).*experiment_data(i).signal_ratio(fit_ind));
    experiment_data(i).scale_factor = 1./experiment_data(i).signal_ratio_times_distance;
end

%% set up model

% spcify constants
g = 9.8;
c = 1480;
rho = 1000;
kappa = 7/5;
sigma = 0.072;
p_atm = 101.3e3;

% tank reflection coefficients
beta_wall = [-0.9 -0.8 -0.7 -0.6];
beta_surface = -1;

% source distances from receiver
distance = [0.0016 0.0025 0.004 0.0063 0.01 0.016 0.025 0.04 0.063 0.1 0.16 0.25 0.35 0.4];

% bubble radius
R_eq = [1.65 1.18]*1e-3; % m
forcing_amplitude = [1000 1000]; % Pa

% define tank size
Lx = 0.552;
Ly = 0.43;
Lz = 0.405;

% specify source location
x_source = 0.15;
y_source = 0.26;
z_source = 0.285;
r_source = [x_source; y_source; z_source];

% compute receiver position
r_receiver = zeros(3,length(distance));
for j = 1:length(distance)
    r_receiver(:,j) = r_source + [distance(j); 0; 0];
end

depth = Lz - z_source;
wall_distance = Lx - x_source;
p_inf = p_atm + rho*g*depth;

% integration properties
integration_time_step = 1e-6;
duration = 10e-3;
t = (0:integration_time_step:duration)';

% compute bubble natural frequency
[natural_frequency,~] = compute_bubble_natural_frequency(R_eq,p_inf,kappa,sigma,rho);

% create synthetic pulses
t0 = 0.2e-3; % s
tau = [5 2.9]*1e-3; % s
f = natural_frequency;
p_radiated_sythetic = sin(2*pi*f.*(t-t0)).*exp(-(t-t0)./tau);
p_radiated_sythetic(t<t0,:) = 0;

%% compute image properties for source and receiver locations
[image_distances_at_source,image_coefficients_at_source] = compute_source_image_distances_and_reflection_coefficients(r_source,r_source,Lx,Ly,Lz,c,beta_wall,beta_surface,duration);
image_distances_at_receiver = zeros(length(image_distances_at_source),1,length(distance));
image_coefficients_at_receiver = zeros(length(image_distances_at_source),length(beta_wall),length(distance));

for k = 1:length(distance)
    [image_distances_at_receiver(:,:,k),image_coefficients_at_receiver(:,:,k)] = compute_source_image_distances_and_reflection_coefficients(r_source,r_receiver(:,k),Lx,Ly,Lz,c,beta_wall,beta_surface,duration);
end

%% initialize arrays
p_radiated_free_bubble = zeros(length(t),length(R_eq));
p_receiver_free_field_bubble = zeros(length(t),length(R_eq),length(distance));
p_receiver_free_field_synthetic = zeros(length(t),length(R_eq),length(distance));

p_receiver_uncoupled = zeros(length(t),length(R_eq),length(beta_wall),length(distance));

p_radiated_coupled = zeros(length(t),length(R_eq),length(beta_wall));
p_receiver_coupled = zeros(length(t),length(R_eq),length(beta_wall),length(distance));

%% free field propagation
% compute radiated signal using RP equation when bubble is not influenced
% by reflections from tank walls
parfor i = 1:length(R_eq)
    [~, ~, ~, ~, p_radiated_free_bubble(:,i), ~] = integrate_rayleigh_plesset_equation_free(integration_time_step,duration,R_eq(i),forcing_amplitude(i),kappa,rho,depth);
end

% compute pressure field at receiver for free field propagation
parfor i = 1:length(distance)
    p_receiver_free_field_bubble(:,:,i) = compute_tank_reflection_time_domain(t,p_radiated_free_bubble,c,distance(i),1);
    p_receiver_free_field_synthetic(:,:,i)  = compute_tank_reflection_time_domain(t,p_radiated_sythetic,c,distance(i),1);
end

%% compute radiated signal when tank reflections influence bubble behaviour
for i = 1:length(R_eq)
    parfor j = 1:length(beta_wall)
        [~, ~, ~, ~, p_radiated_coupled(:,i,j), ~] = integrate_rayleigh_plesset_equation_coupled(integration_time_step,duration,R_eq(i),forcing_amplitude(i),kappa,rho,c,depth,image_distances_at_source,image_coefficients_at_source(:,j));
    end
end

%% compute received signal with tank reflection
for j = 1:length(beta_wall)
    parfor k = 1:length(distance)
        p_receiver_uncoupled(:,:,j,k) = compute_tank_reflection_time_domain(t,p_radiated_sythetic,c,image_distances_at_receiver(:,:,k),image_coefficients_at_receiver(:,j,k));
        p_receiver_coupled(:,:,j,k) = compute_tank_reflection_time_domain(t,p_radiated_coupled(:,:,j),c,image_distances_at_receiver(:,:,k),image_coefficients_at_receiver(:,j,k));
    end
end

%% compute statistics
rms_amplitude_free_field_bubble = zeros(length(R_eq),length(distance));
rms_amplitude_free_field_synthetic = zeros(length(R_eq),length(distance));

[rms_amplitude_coupled,rms_amplitude_uncoupled] = deal(zeros(length(R_eq),length(beta_wall),length(distance)));

for i = 1:length(R_eq)
    for k = 1:length(distance)
        rms_amplitude_free_field_bubble(i,k) = sqrt(mean(p_receiver_free_field_bubble(:,i,k).^2));
        rms_amplitude_free_field_synthetic(i,k) = sqrt(mean(p_receiver_free_field_synthetic(:,i,k).^2));
    end
end

for i = 1:length(R_eq)
    for j = 1:length(beta_wall)
        for k = 1:length(distance)
            rms_amplitude_uncoupled(i,j,k) = sqrt(mean(p_receiver_uncoupled(:,i,j,k).^2));
            rms_amplitude_coupled(i,j,k) = sqrt(mean(p_receiver_coupled(:,i,j,k).^2));
        end
    end
end

[model_scale_factor_coupled,model_scale_factor_uncoupled] = deal(zeros(length(R_eq),1));
for i = 1:length(R_eq)
    model_scale_factor_coupled(i) = 1./(rms_amplitude_free_field_bubble(i,1).*distance(1));
    model_scale_factor_uncoupled(i) = 1./(rms_amplitude_free_field_synthetic(i,1).*distance(1));
end

%% plot statistics
cmap = cbrewer2('set1',length(beta_wall)+1);
filled_size = 20;
unfilled_size = 40;

for i = 1:length(R_eq)
    figure(i);
    clf;
    hold on;
    h(1) = scatter(NaN,NaN,filled_size,[0 0 0],'filled','displayname','Filled = coupled');
    h(2) = scatter(NaN,NaN,unfilled_size,[0 0 0],'displayname','Open = uncoupled');
    h(length(beta_wall)+3) = scatter(experiment_data(i).distance,experiment_data(i).scale_factor.*experiment_data(i).signal_ratio,unfilled_size,cmap(j+1,:),'s','displayname','Tank data'); %#ok<SAGROW>
    h(length(beta_wall)+4) = plot(distance,model_scale_factor_coupled(i).*rms_amplitude_free_field_bubble(i,:),'color','black','displayname','Free field scaling'); %#ok<SAGROW>
    for j = 1:length(beta_wall)
        scatter(distance,model_scale_factor_uncoupled(i).*squeeze(rms_amplitude_uncoupled(i,j,:)),unfilled_size,cmap(j,:));
    end
    for j = 1:length(beta_wall)
        h(j+2) = scatter(distance,model_scale_factor_coupled(i).*squeeze(rms_amplitude_coupled(i,j,:)),filled_size,cmap(j,:),'filled','displayname',['Model, \beta = ' num2str(beta_wall(j))]);
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlims = get(gca,'xlim');
    set(gca,'xlim',[xlims(1) wall_distance])
    hold off;
    grid on;
    xlabel('Bubble hydrophone distance (m)');
    ylabel('Gated pulse RMS amplitude ratio');
    title(['Natural frequency: ' num2str(natural_frequency(i)/1e3,'%.1f') ' kHz']);
    legend(h,'location','southwest');
    exportgraphics(gcf,['tank_model_vs_experiment_' num2str(f(i),'%.0f') '_hz.pdf']);
end

%% plot time series
%{
figure(3);
plot(t,p_radiated_free(:,1));
hold on;
plot(t,p_radiated_coupled(:,1,1))
hold off;

figure(4)
plot(t,p_receiver_uncoupled(:,1,1,1))

%}