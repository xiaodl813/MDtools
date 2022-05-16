clear; close all; clc;

%% Nanotube parameters

%%%%%%%%%%%%%%% CHANGE HERE %%%%%%%%%%%%%%%

Lz_value = 120; % unit is layer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% chiral index (n,m)
Chiral_n = 10;
Chiral_m = 10;

% Constant of nanotube
Bond_Length = 1.432; % average bond length, unit is Angstrom
Layer_Length = 2.475; % unit is Angstrom (per layer)

% calulate volume of nanotube
a = sqrt(3)*Bond_Length; % unit vector of nanotube, unit is Angstrom
C = a * sqrt(Chiral_m^2 + Chiral_n^2 + ...
    Chiral_n * Chiral_m); % perimeter of nanotube, unit is Angstrom
Lz = Lz_value * Layer_Length; % length of the whole nanotube, unit is A

atom_r = Bond_Length /2; % average atomic radius, unit is Angstrom
% NOTE: `atom_r` is assumed to be half the average bond length

V = C * Lz * (2 * atom_r); % unit is A^3

%% Read file

lp = [1 2 3 4 5 10 20 30 60];

for i = lp
    eval(['file_' num2str(i) '= "hac' num2str(i) '.txt"']);
    eval(['[hac_' num2str(i) ', parameter] = read_hac(file_' num2str(i) ');']);
    % NOTE: for every i, `parameter_i` is equal
end

%% Compute kappa (Green-Kubo)

% parameters of Green-Kubo formula
dt = 0.0005; % unit is ps
T = 300; % unit is K
kB = 8.617e-5; % unit is eV/K

% calculate kappa
scale = 1.6e3 / (kB * T^2 * V) * (parameter.Ns * dt);
t = (1 : parameter.Nc) * parameter.Ns * dt; % correlation time

RTC = zeros(20,9); % kappa at t=200ps
Kappa = zeros(2,9); % Kappa = [Kappa_i, Kappa_err_i]

for i = 1 : size(lp, 2)
    eval(['rtc_' num2str(lp(i)) ...
          '= scale * cumtrapz(hac_' num2str(lp(i)) ');']);
    RTC(:, i) = eval(['rtc_' num2str(lp(i)) '(end, :)'';']);

    % suppose that kappa converges after 100ps
    eval(['kappa_con_' num2str(lp(i)) ' = mean(rtc_' num2str(lp(i)) ...
          '(parameter.Nc/2 + 1 : parameter.Nc, :));']);
    Kappa(1, i) = eval(['mean(kappa_con_' num2str(lp(i)) ');']); % Kappa_i
    Kappa(2, i) = eval(['std(kappa_con_' num2str(lp(i)) ...
                        ') / sqrt(parameter.M);']); % Kappa_err_i
end
% t-HFACF & t-kappa

for i = lp
    % t-HFACF
    eval(['plot_HFACF(t,hac_' num2str(i) ',"HFACF_' num2str(i) '.jpg");']);
    % t-kappa
    eval(['plot_kappa(t,rtc_' num2str(i) ',"kappa_' num2str(i) '.jpg");']);
end

%% boxplot
% NOTE: the plot code should be packaged

font_size = 15;

figure
boxplot(RTC,'Notch','on','Labels',{1,2,3,4,5,10,20,30,60});
hold on
scatter(1:9,Kappa(1,:));
xlabel('$\ell_p$ (layers)', 'fontsize', font_size, 'interpreter','latex');
ylabel('$\kappa\ ({\rm Wm^{-1}K^{-1}})$', 'fontsize', font_size,...
                                          'interpreter','latex');
set(gca, 'fontsize', font_size);
saveas(gcf, "boxplot.jpg");
close

%% scatter plot (accurater kappa fit by t-distribution)
% NOTE: the plot code should be packaged

Kappa_fit = zeros(3,9);
for i = 1 : size(lp,2)
    eval(['nd_' num2str(lp(i)) ' = fitdist(RTC(:,i), ''Normal'');']);
    eval(['Kappa_fit(1, i) = nd_' num2str(lp(i)) '.mu;']);
    format shortG
    eval(['Kappa_95CI = paramci(nd_' num2str(lp(i)) ');']);
    Kappa_fit(2:3, i) = Kappa_95CI(:, 1);
end

figure
scatter(lp,Kappa_fit(1,:),'x');
hold on
plot([1,2,3,5,10,20,60],Kappa_fit(1,[1:3 5:7 9]),'--');
errorbar(lp,Kappa_fit(1,:),Kappa_fit(1,:)-Kappa_fit(2,:),Kappa_fit(3,:)-Kappa_fit(1,:),'b');
xlabel('$\ell_p$ (layers)', 'fontsize', font_size, 'interpreter','latex');
ylabel('$\kappa\ ({\rm Wm^{-1}K^{-1}})$', 'fontsize', font_size,...
                                          'interpreter','latex');
set(gca, 'fontsize', font_size);

%%

