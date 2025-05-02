classdef Parameters
    % Define constant variables used throughout the global function
    properties (Constant = true)  
    %------------------ Physical constants in cgs units------------------------
    c=299792458;                                    %ms^-1;
    G=6.6740831*10^-11  ;                           %m^3*kg-1
    Ms=1.988489*10^30 ;                          %kg
    pc=3.08567758149*10^16 ;                        %m
    AU=149597870700;
    H0 = 72/(3.*10.^19) ;                      % Hubble parameters in sec, cgs 
    Gpc = 1.e9*3.0856776*1e16;                 % GParsec, cgs 
    T1year = 1.*3.127.*10.^7;
    R      = 499;
 %------------------------------------------------------------------------   
    end
    methods

        function MK = readinput(obj,filename)
            MK = readmatrix(filename, 'Delimiter', ','); 
        end
        
        function ret=distance(obj,z)
                %      This is the EM luminosity distance.
                % 
                % Parameters
                % ----------
                % z : scalar, redshift to source
                % Omegam0, OmegaLAMB0, H0: Cosmological parameters today (z = 0)
                % 
                % Returns
                % -------
                % dL: scalar, luminosity distance to the source
                    O_M=  0.29;
                    O_A = 0.71;     
                    zz=linspace(0,z,50);dzz=zz(3)-zz(2);  
                    ret=obj.c./obj.H0.*(1+z).*sum(1./sqrt(O_M.*(1+zz).^3+O_A)).*dzz; %Bennet C L et al 2003 Astrophys. J. S. 148 1
        end 
    end
    
end

