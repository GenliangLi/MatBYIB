classdef Detector
    %DETECTOR  Class implementing gravitational wave detector models
    %   Provides detector response functions and noise models for:
    %     - LISA (Laser Interferometer Space Antenna)
    %     - LIGO (Laser Interferometer Gravitational-Wave Observatory) 
    %     - Taiji (Chinese space-based detector)
    %
    %  Methods:
    %    Response_function    - Compute antenna pattern responses
    %    PowerSpectralDensity - Generic PSD model
    %    Sen_curve_LISA       - LISA sensitivity curve
    %    Sen_curve_LIGO       - LIGO sensitivity curve  
    %    Sen_curve_TAIJI      - Taiji sensitivity curve
    
    properties (Constant = true)  
        % Note: Add any detector constants here (e.g., arm lengths)
    end
    
    methods       
        function [F1plus,F1cross,F2plus,F2cross] = Response_function(obj,theta,phi,psi)
            %RESPONSE_FUNCTION  Compute detector antenna pattern responses
            %
            %  Inputs:
            %    theta : Source polar angle in ecliptic coordinates [rad]
            %    phi   : Source azimuthal angle in ecliptic coordinates [rad] 
            %    psi   : Wave polarization angle [rad]
            %
            %  Outputs:
            %    F1plus, F1cross : Response factors for + and × polarizations (Channel 1)
            %    F2plus, F2cross : Response factors for + and × polarizations (Channel 2)
            %
            %  Equations:
            %    Implements Eqs. (7)-(8) from arXiv:1702.00786 (LISA Science Case)
            
            F1plus  = ((1 + cos(theta).^2).*cos(2.*phi).*cos(2.*psi))/2 ...
                     - cos(theta).*sin(2.*phi).*sin(2.*psi);
            F1cross = ((1 + cos(theta).^2).*cos(2.*phi).*sin(2.*psi))/2 ...
                     + cos(theta).*sin(2.*phi).*cos(2.*psi);

            F2plus  = ((1 + cos(theta).^2).*sin(2.*phi).*cos(2.*psi))/2 ...
                     + cos(theta).*cos(2.*phi).*sin(2.*psi);
            F2cross = ((1 + cos(theta).^2).*sin(2.*phi).*sin(2.*psi))/2 ...
                     - cos(theta).*cos(2.*phi).*cos(2.*psi);
        end
        
        function [Sn] = PowerSpectralDensity(obj, f)
            %POWERSPECTRALDENSITY  Generic PSD model for validation
            %
            %  Inputs:
            %    f : Frequency vector [Hz]
            %
            %  Outputs:
            %    Sn : Power spectral density [Hz^(-1)]
            %
            %  Notes:
            %    Reference model from Cornish & Robson (2017)
            %    DOI: 10.1088/1361-6382/aa5f50
            
            f0 = 215;       % Reference frequency [Hz]
            S0 = 1e-49;     % Reference amplitude [Hz^(-1)]
            x = f ./ f0;    % Normalized frequency
            
            Sn = S0 .* (x.^(-4.14) - 5.*x.^(-2) + ...
                (111.*(1 - x.^2 + x.^4/2))./(1 + x.^2/2));
        end
        
        function [Sn] = Sen_curve_LISA(obj, f)
            %SEN_CURVE_LISA  Compute LISA sensitivity curve
            %
            %  Inputs:
            %    f : Frequency vector [Hz]
            %
            %  Outputs:
            %    Sn : Strain sensitivity [Hz^(-1/2)]
            %
            %  Implementation:
            %    Combines acceleration noise (P_acc), optical metrology noise (P_OMS),
            %    and galactic confusion noise (Sc) following:
            %    - Eqs. (1)-(3) from LISA Science Requirements Document
            %    - Arm length L = 2.5e9 m
            
            c = 299792458;   % Speed of light [m/s]
            L = 2.5e9;       % Arm length [m]
            f_star = c/(2*pi*L);  % Transfer frequency [Hz]
            
            % Acceleration noise (proof mass noise)
            P_acc = (3e-15)^2 .* (1 + (0.4e-3./f).^2);  % [m^2 s^(-4) Hz^(-1)]
            
            % Optical Metrology System noise
            P_OMS = (1.5e-11)^2;  % [m^2 Hz^(-1)]
            
            % Total noise components
            Pn = P_OMS/L^2 + 2*(1 + cos(f/f_star).^2).*P_acc./(2*pi*f).^4/L^2;
            
            % Galactic confusion noise parameters
            as = 0.138;      % Exponential decay coefficient
            betas = -221;    % Sinusoidal modulation coefficient  
            keta = 521;      % Modulation frequency coefficient
            gamma = 1680;    % Transition width coefficient
            f_k = 0.00113;   % Transition frequency [Hz]
            A = 9e-45;       % Amplitude [Hz^(-7/3)]
            
            % Galactic binary confusion noise
            Sc = A*f.^(-7/3).*exp(-f.^as + betas*f.*sin(keta*f)).*...
                (1 + tanh(gamma*(f_k - f)));
            
            % Combined sensitivity curve
            Sn = 10/(3*L^2)*(P_OMS + 2*(1 + cos(f/f_star).^2).*...
                P_acc./(2*pi*f).^4).*(1 + 0.6*(f/f_star).^2) + Sc;
            
            % Ensure smooth output via spline interpolation
            Sn = interp1(f, Sn, f, 'spline');  
        end
        
        function [Sn] = Sen_curve_LIGO(obj, F)
            %SEN_CURVE_LIGO  Compute LIGO design sensitivity
            %
            %  Inputs:
            %    F : Frequency vector [Hz]
            %
            %  Outputs:
            %    Sn : Strain sensitivity [Hz^(-1/2)]
            %
            %  Model:
            %    Analytic approximation to LIGO O4 design sensitivity
            %    Reference: LIGO-T1800044-v5
            
            f0 = 70;         % Knee frequency [Hz]
            S0 = 3e-48;      % Normalization [Hz^(-1)]
            Sn = S0*((f0./F).^4 + 2 + 2*(F./f0).^2);
        end

        function [Snfn] = Sen_curve_TAIJI(obj, f)
            %SEN_CURVE_TAIJI  Compute Taiji sensitivity from data file
            %
            %  Inputs:
            %    f : Frequency vector [Hz]
            %
            %  Outputs:
            %    Snfn : Interpolated strain sensitivity [Hz^(-1/2)]
            %
            %  Data Source:
            %    TAIJISEN.txt - Official sensitivity curve from Taiji Mission
            %    Reference: arXiv:1807.09495
            
            sp = load('TAIJISEN.txt');  % Load [freq(Hz), sqrt(PSD)]
            ff = sp(:,1);               % Frequency grid
            Snf = sp(:,2).^2;           % Convert to PSD
            Snfn = interp1(ff, Snf, f, 'spline');  % Cubic spline interpolation
        end 
    end
end
    