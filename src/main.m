clear
clc
%%%%%% Parameter values for fiducial model. 
% Add/remove parameters as needed, but ensure all other relevant parts of the code are adjusted. 

% Load physical constants and unit conversions from Common_CF module
Units = Common_CF;

% Read input parameters from 'input.txt' file
parameters = Units.readinput('input.txt');

% Extract and assign physical parameters with units and descriptions:
M0    = parameters(1);     % Central black hole mass [solar masses]
mu0   = parameters(2);     % Orbiting object mass [solar masses]
e_LSO = parameters(3);     % Eccentricity at last stable orbit (dimensionless)
spin0 = parameters(4);     % Dimensionless spin parameter of central BH [S/M^2]
Zf    = parameters(5);     % Source redshift (affects luminosity distance)

% Calculate luminosity distance using cosmological relations
DL = Units.distance(Zf);

%----------------------------- Geometrical parameters and initial conditions ------------------------------
% Convert all angular parameters from degrees to radians:
lambda = parameters(6).*pi/180;     % Angle between orbital and spin angular momentum vectors
Phi_LSO = parameters(7).*pi/180;    % Mean anomaly at last stable orbit
Gamma_LSO = parameters(8).*pi/180;  % Angle between L×S and pericenter
alpha_LSO = parameters(9).*pi/180;  % Azimuthal direction of L in orbital plane
thetaS = parameters(10).*pi/180;    % Source polar angle in ecliptic coordinates
phiS = parameters(11).*pi/180;      % Source azimuthal angle in ecliptic coordinates
thetaK = parameters(12).*pi/180;    % Spin vector polar angle
phiK = parameters(13).*pi/180;       % Spin vector azimuthal angle
phi0 = parameters(14).*pi/180;      % Initial orbital phase
t_max = parameters(15);              % Maximum evolution time [seconds]

%--------------------------------------------------------------------------------------------------------
% Package all parameters into a single array for waveform generation
params = [M0,mu0,e_LSO,spin0,Zf,lambda,Phi_LSO,Gamma_LSO,alpha_LSO,thetaS,phiS,thetaK,phiK,phi0,t_max];

% Load AK waveform module and Initialize AK waveform generator with the parameters
AKwave = AKwaveform(params);

% Calculate orbital evolution (returns time series of orbital parameters)
[t, Phi, nu, e, gamma, alpha] = AKwave.evalution;

% Load plotFigure module and plot orbital parameter evolution
plotFigure = plotFigure;
plotFigure.plotOrbit(t, nu, e, Phi);

% Generate time-domain gravitational waveform
[t, ht] = AKwave.waveform_td;

% Compute Fourier transform to get frequency-domain waveform
[fsg, hf] = AKwave.Fourier_tran(t, ht);

% Plot time-domain waveform
plotFigure = plotFigure;
plotFigure.plot_t_fig(t, ht, 0);

% Plot frequency-domain waveform
plotFigure.plot_f_fig(fsg, hf, 0);

% Create finer frequency bin array for interpolation
freq_bin = linspace(0.1*nu(1), nu(end)*11, 1e6); % Logarithmic frequency bins

% Interpolate waveform onto new frequency bins
hff = interp1(fsg, hf, freq_bin, 'nearest');
% Replace any NaN values with zeros (edge cases)
hff(find(isnan(hff))) = 0;

% Plot interpolated frequency-domain waveform ,
%"0" repreasents a logarithmic plot of the frequency-domain waveform
%"1" repreasents a nature plot of the frequency-domain waveform
plotFigure.plot_f_fig(freq_bin, hff, 0);


% Load Fisher Matrix module
FM = Fisher_Matrix;    

% Calculate frequency bin width (uniform spacing assumed)
delta_f = freq_bin(2) - freq_bin(1);

% Load Fisher detector module object and get LISA sensitivity curve
Det = Detector;         
PSD = Det.Sen_curve_LISA(freq_bin); 

% Set number of parameters for Fisher Matrix calculation (NN=4 for basic params)
NN =2;  % Typically mass1, mass2, eccentricity, spin

% Calculate covariance matrix using Fisher Matrix formalism
[cov_Matrix] = FM.Matrix(params, freq_bin, PSD, NN);

% Interpolate waveform onto frequency bins (clean NaN values)
h_f = interp1(fsg, hf, freq_bin, 'nearest');
h_f(isnan(h_f)) = 0;

% Calculate and display signal-to-noise ratio (SNR)
SNR = FM.inner_prod(h_f, h_f, PSD, delta_f);  % SNR squared calculation
str = ['SNR for source: ', num2str(sqrt(SNR))]; % Convert to linear SNR
disp(str);  % Display SNR in command window

% %---------------------------------------------------------
% % Fisher Matrix Calculation Breakdown (alternative method)
% %---------------------------------------------------------
% 
% % Calculate parameter derivatives numerically
% for ii = 1:NN
%     % Numerical derivative of waveform with respect to parameter ii
%     diff_vec(ii,:) = FM.diff_param(ii, params, freq_bin); 
% end
% 
% 
% % Clean numerical derivatives (remove NaN values)
% diff_vec(isnan(diff_vec)) = 0;
% 
% 
% % Initialize Fisher Matrix (identity matrix as placeholder)
% fish_mix = eye(NN);  % NN x NN identity matrix
% 
% % Build Fisher Matrix through inner products of parameter derivatives and
% % calculate FM elements using noise-weighted inner product
% for ii = 1:NN
%     for jj = 1:NN
%         fish_mix(ii,jj) = FM.inner_prod(diff_vec(ii,:), diff_vec(jj,:), PSD, delta_f);
%     end
% end
% 
% % Calculate covariance matrix as inverse of Fisher Matrix
% fish_mix_inv = inv(fish_mix); 
% 
% % Store covariance matrix (alternative to direct output from FM.Matrix)
% cov_Matrix = eye(NN);  % Initialize
% for ii = 1:NN
%     for j = 1:NN
%         cov_Matrix(ii,j) = fish_mix_inv(ii,j);  % Copy inverted matrix
%     end
% end
% 
% Display parameter uncertainties (1-sigma errors)
for ii = 1:NN
    dparams(ii) = sqrt(cov_Matrix(ii,ii));  % Standard deviation = sqrt(variance)
    str = ['Parameter', num2str(ii), ' ','uncertainty: ', num2str(dparams(ii))];
    disp(str);  % Print each parameter's uncertainty
end

%---------------------------------------------------------
% MCMC (Markov Chain Monte Carlo) Implementation for Gravitational Wave Parameter Estimation
%---------------------------------------------------------
% Calculate noise variances for each frequency bin
% PSD is the power spectral density, delta_f is frequency bin width
% Variance = PSD/(4Δf) for likelihood calculation
variances = (PSD)./(4.*delta_f); 

% Generate random noise realizations (real and imaginary parts)
noise = normrnd(0, sqrt(variances));      % Real component of noise
noise_j = normrnd(0, sqrt(variances));    % Imaginary component of noise

% Construct complex frequency-domain noise
noise_f = noise + 1j.*noise_j; 

% Visualize noise and signal for debugging
figure; plot(freq_bin, abs(noise_f), '*'); % Plot noise realization
figure; plot(freq_bin, abs(hff));          % Plot signal template

% Create simulated data
data_f = h_f + noise_f;

% Start timer for performance monitoring
tic;

% MCMC Configuration
Ntotal = 2000;    % Total number of MCMC iterations
burnin = 0;       % Burn-in period (0 since we start at true parameters)

%--------------------------------------------------------------------------
% Define Uniform Prior Distributions for Parameters
%--------------------------------------------------------------------------

% Central black hole mass [M_sun]
M0_low = 1.e+5;  
M0_high = 1.e+7;    

% Orbiting object mass [M_sun]
m0_low = 1.e+0;
m0_high = 1.e+2;

% Black hole spin parameter [dimensionless]
spin0_low = 0.001;
spin0_high = 0.1;

% Redshift
Zf_low = 0.001;
Zf_high = 1;

% Eccentricity at last stable orbit
e_LSO_low = 0.2;
e_LSO_high = 0.4;

% Angular parameters [radians]
lambda_low = 0;
lambda_high = pi/2;

Phi_LSO_low = 0;
Phi_LSO_high = pi/2;

Gamma_LSO_low = 0;
Gamma_LSO_high = pi/2;

alpha_LSO_low = 0;
alpha_LSO_high = pi/2;

thetaS_low = 0;
thetaS_high = pi/2;

% Separate parameters into those being sampled (first NN) and fixed parameters
true_vals = params(1:NN);  % Parameters to estimate (first NN parameters)
params_sh = params(NN+1:end); % Fixed parameters (held constant)

%--------------------------------------------------------------------------
% Run MCMC Sampling
%--------------------------------------------------------------------------
printerval = 200;  % Interval for printing progress updates

% Load MCMC module
MCMC = MCMC; 

% Run main MCMC sampling routine
[params_samples, lp] = MCMC.MCMC_run(data_f, freq_bin, true_vals, cov_Matrix,...
                                     variances, params_low, params_sh, params_high,...
                                     Ntotal, printerval);
% Display elapsed time
toc;

%--------------------------------------------------------------------------
% Convergence Diagnostics and Continued Sampling
%--------------------------------------------------------------------------
% Check convergence using Gelman-Rubin statistic (R)
[samples_tol, lps_tol, R] = MCMC.converg(params_samples, lp);

% Continue extending chains until all R < 1.05 (convergence criterion), 
while any(R(:) > 1.05)
    Ntotal = floor(1.2*Ntotal);
    [sam_con, lp_con] = MCMC.MCMC_contin(data_f, freq_bin, cov_Matrix, variances,...
                                        params_low, params_sh, params_high,...
                                        sam_con, lp_con, Ntotal,...
                                        accept_reject_count, printerval);
    [samples_tol, lps_tol, R] = MCMC.converg(sam_con, lp_con);
end

% Set up the figure
samples_MC = samples_tol;
Names ={'m_1','m_2','e_{LSO}','S'};
labels ={'m_1','m_2','e_{LSO}','S'};

plotFigure=plotFigure;
% Generate samples
samples_FM= mvnrnd(true_vals, cov_Matrix, 50000);

for i = 1:length(Names)
    median_val = median(samples_FM(:, i));
    lower_val = median_val - prctile(samples_FM(:, i), 16);
    upper_val = prctile(samples_FM(:, i), 84) - median_val;
    titles{i} = plotFigure.format_title(Names{i}, median_val, lower_val, upper_val);
end
fig = plotFigure.cornerplot(samples_FM(:,1:NN),labels,true_vals,titles);

for i = 1:length(Names)
    median_val = median(samples_MC(:, i));
    lower_val = median_val - prctile(samples_MC(:, i), 16);
    upper_val = prctile(samples_MC(:, i), 84) - median_val;
    titles{i} = plotFigure.format_title(Names{i}, median_val, lower_val, upper_val);
end
fig = plotFigure.cornerplot(samples_MC(:,1:NN),labels,true_vals,titles);

str=['Source'];
disp(str);
str=['MCMC M0:',num2str(sqrt(var(samples_MC(:,1))))];
disp(str);
str=['MCMC mu0:',num2str(sqrt(var(samples_MC(:,2))))];
disp(str);
str=['MCMC e_LSO:',num2str(sqrt(var(samples_MC(:,3))))];
disp(str);
str=['MCMC spin0:',num2str(sqrt(var(samples_MC(:,4))))];
disp(str);

sdh=sqrt(diag(cov_Matrix));
str=['FM:',num2str(sdh(1)),' ',num2str(sdh(2)),' ',num2str(sdh(3)),' ', num2str(sdh(4))];
disp(str);



