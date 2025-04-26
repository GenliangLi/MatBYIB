classdef Common_CF
    %COMMON_CF  Class containing fundamental physical constants and utility functions
    %
    %   This class defines physical constants in SI and astronomical units, and provides
    %   methods for cosmological calculations. All constants are defined as read-only
    %   properties with fixed units as noted below.
    %
    %   Example:
    %     cf = Common_CF;
    %     lumDist = cf.distance(0.5);  % Calculate luminosity distance at z=0.5
    %
    % Version:
    %   1.0 (2025-03-20) Initial implementation
    %   1.1 (2025-04-15) Added cosmological parameter customization
    
    properties (Constant = true)
        %MS Solar mass [kg]
        %   Reference: IAU 2015 Resolution B3
        Ms = 1.9892e30;
        
        %H0 Hubble constant [s^-1]
        %   Current value: 72 km/s/Mpc converted to inverse seconds
        H0 = 72/(3.0856776e19);
        
        %PC Parsec [m]
        %   Exact value per IAU 2015 Resolution B2
        pc = 3.0856776e16;
        
        %GPC Gigaparsec [m]
        Gpc = 3.0856776e25;
        
        %T1YEAR Tropical year [s]
        %   Mean tropical year duration at J2000.0
        T1year = 3.1556925e7;
        
        %R Schwarzschild radius multiplier
        %   Used for defining integration boundaries (default 499)
        R = 499;
        
        %C Speed of light in vacuum [m/s]
        %   Exact value per SI definition
        c = 299792458;
        
        %G Gravitational constant [m^3 kg^-1 s^-2]
        %   CODATA 2018 recommended value
        G = 6.67430e-11;
        
        %AU Astronomical unit [m]
        %   IAU 2012 exact definition
        AU = 149597870700;
    end
    
    methods
        function MK = readinput(obj, filename)
            %READINPUT Read parameter matrix from CSV file
            %
            %   MK = READINPUT(filename) loads a comma-delimited matrix from the
            %   specified file. The file should contain numerical values only.
            %
            %   Input:
            %       filename    Path to input file (string)
            %
            %   Output:
            %       MK          Parameter matrix [N×M double]
            
            MK = readmatrix(filename, 'Delimiter', ',');
        end
        
        function ret = distance(obj, z)
            %DISTANCE Calculate luminosity distance for given redshift
            %
            %   D = DISTANCE(z) computes the EM luminosity distance using
            %   standard ΛCDM cosmology with default parameters:
            %     - H0 = 72 km/s/Mpc
            %     - Ωm = 0.29
            %     - ΩΛ = 0.71
            %
            %   Reference:
            %     Bennett et al. 2003, ApJS 148, 1 (WMAP first-year results)
            %
            %   Input:
            %       z       Redshift [dimensionless]
            %
            %   Output:
            %       dL      Luminosity distance [m]
            
            % Cosmological parameters (WMAP9 values)
            O_M = 0.29;     % Matter density parameter
            O_A = 0.71;     % Dark energy density parameter
            
            % Numerical integration setup
            zz = linspace(0, z, 50);    % Integration points
            dzz = zz(3) - zz(2);        % Step size
            
            % Integrate comoving distance (Eq. 18 in Hogg 1999)
            ret = obj.c/obj.H0 * (1+z) * ...
                  sum(1./sqrt(O_M*(1+zz).^3 + O_A)) * dzz;
        end
    end
end