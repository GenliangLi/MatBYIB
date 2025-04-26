 classdef waveform
   % write a description of the class here.
       properties 
       % define the properties of the class here, (like fields of a struct)
            M0;   
            mu0;
            e_LSO ;
            spin0;
            Zf; 
            lambda ;
            Phi_LSO ;
            Gamma_LSO ;
            alpha_LSO ;
            thetaS;
            phiS;
            thetaK;
            phiK;
            phi0;
            t_max;
            n_max   = 11  ;                                    % maximum number of overtones to include in the waveform. 
       end
      
       methods
       % methods, including the constructor are defined in this block
           
           function obj = waveform(parmas)
           % class constructor
               if(nargin > 0)
                obj.M0           =   (parmas(1)) ;  
                obj.mu0          =   parmas(2) ; 
                obj.e_LSO        =   parmas(3) ;
                obj.spin0        =   parmas(4) ; 
                obj.Zf           =   parmas(5) ; 
                obj.lambda       =   parmas(6) ; 
                obj.Phi_LSO       =   parmas(7) ; 
                obj.Gamma_LSO     =   parmas(8) ; 
                obj.alpha_LSO     =   parmas(9) ; 
                obj.thetaS       =   parmas(10) ; 
                obj.phiS         =   parmas(11) ; 
                obj.thetaK       =   parmas(12) ; 
                obj.phiK         =   parmas(13) ; 
                obj.phi0         =   parmas(14) ; 
                obj.t_max        =   parmas(15) ; 
               end
           end
           
           function [hff]=waveform_fd(obj,freq_bin) 
           
                    [t,ht]=waveform_td(obj);
                  
                    [fsg,hf]   =Fourier_tran(obj,t,ht);

                    hff     = interp1(fsg,hf,freq_bin,'nearest');
           end   
           
           function [t,ht]=waveform_td(obj)
                
                [t,Phi,nu,e,gamma,alpha]=evalution(obj);
                
                x = ((2.*pi.*Common_CF.G.*obj.M0*Common_CF.Ms.*nu)/Common_CF.c.^3);
                A = (x.^(2/3).*(Common_CF.G.*obj.mu0.*Common_CF.Ms/Common_CF.c.^2));
                thetaL = acos(cos(obj.thetaK).*cos(obj.lambda) + sin(obj.thetaK).*sin(obj.lambda).*cos(alpha)) ;
                SB1    = (sin(obj.thetaK).*sin(obj.phiK).*cos(obj.lambda) - cos(obj.phiK).*sin(obj.lambda).*sin(alpha) - sin(obj.phiK).*cos(obj.thetaK).*sin(obj.lambda).*cos(alpha));
                SB2    = (sin(obj.thetaK).*cos(obj.phiK).*cos(obj.lambda) + sin(obj.phiK).*sin(obj.lambda).*sin(alpha) - cos(obj.phiK).*cos(obj.thetaK).*sin(obj.lambda).*cos(alpha));
                phiL    = atan2(real(SB1),real(SB2));

                Ln = cos(obj.thetaS).*cos(thetaL) + sin(obj.thetaS).*sin(thetaL).*cos(obj.phiS - phiL);
                Sn = cos(obj.thetaS).*cos(obj.thetaK) + sin(obj.thetaK).*sin(obj.thetaS).*cos(obj.phiS - obj.phiK);
                SB1= (Ln.*cos(obj.lambda) - Sn);
                SB2 = (sin(obj.thetaS).*sin(obj.phiS - obj.phiK).*sin(obj.lambda).*cos(alpha) + (Sn.*cos(obj.thetaK) - cos(obj.thetaS))/(sin(obj.thetaK)).*sin(obj.lambda).*sin(alpha));
                beta = atan2(real(SB1),real(SB2));
                gammatilde = gamma + beta;

                theta = acos(cos(obj.thetaS)/2 - (sqrt(3).*sin(obj.thetaS)/2.*cos((2.*pi.*t)/Common_CF.T1year - obj.phiS))) ;
                SB1   = (sqrt(3).*cos(obj.thetaS) + sin(obj.thetaS).*cos((2.*pi.*t)/Common_CF.T1year - obj.phiS));
                SB2   = (2.*(sin(obj.thetaS).*(sin((2.*pi.*t)/Common_CF.T1year - obj.phiS))));
                phi = (2.*pi.*t)/Common_CF.T1year + atan2(real(SB1), real(SB2));
                SB1 = (cos(thetaL)/2 - (sqrt(3).*sin(thetaL).*cos((2.*pi.*t)/Common_CF.T1year - phiL))/2 - cos(theta).*(cos(obj.thetaS).*cos(thetaL) + cos(phiL - obj.phiS).*sin(obj.thetaS).*sin(thetaL)));
                SB2 = (0.5.*(sin(thetaL).*sin(obj.thetaS).*sin(phiL - obj.phiS)) - 0.5.*(sqrt(3).*cos((2.*pi.*t)/Common_CF.T1year).*(cos(thetaL).*sin(obj.thetaS).*sin(obj.phiS) - cos(obj.thetaS).*sin(thetaL).*sin(phiL))) - 0.5.*(sqrt(3).*sin((2.*pi.*t)/Common_CF.T1year).*(cos(obj.thetaS).*cos(phiL).*sin(thetaL) - cos(obj.phiS).*cos(thetaL).*sin(obj.thetaS))));
                psi = atan2(real(SB1),real(SB2));
                
                Det = Detector;
                [F1plus, F1cross, F2plus, F2cross] = Det.Response_function(theta,phi,psi); % Detector response functions.   
                 
                z = e;
                A_plus_total=0;
                A_cross_total=0;
                for n=1:obj.n_max% sum over orbital overtones n. Integer nmax is the maximum overtone to sum.    
                [a_n, b_n, c_n ]= get_an_bn_cn(obj,A,n,e,Phi);% coefficients entering the waveform.
                [A_plus_n, A_cross_n ]= get_Aplus_A_cross(obj,Ln,a_n,b_n,c_n,gammatilde);% A_plus, A_cross
                A_plus_total =A_plus_total+ A_plus_n;% sum over eigenfrequencies up to nmax for A_plus.
                A_cross_total =A_cross_total+ A_cross_n;% sum over eigenfrequencies up to nmax for A_cross.
                end
                hI  = (0.5.*sqrt(3)).*(F1plus.*A_plus_total + F1cross.*A_cross_total);% waveform with detector response function I
                hII = (0.5.*sqrt(3)).*(F2plus.*A_plus_total + F2cross.*A_cross_total);% waveform with detector response function II

                % total waveform with detector response function without the factor 1/distance in front of it. The prefactor is added conventionaly in the likelihood; the reason being that the prefactor may have parameters to be varied in the McMc. See function iterate_mcmc() below.            
                ht   = (hI + hII)./distance(obj,obj.Zf);  
           end
          
           function [t,Phi,nu,e,gamma,alpha]=evalution(obj)
                   % nu_LSO: The initial condition of the orbital frequency (?) at the LSO. It depends on the central mass and eccentricity.
                    nu_LSO = (Common_CF.c.^3/(2.*pi.*Common_CF.G.*obj.M0.*Common_CF.Ms)).*((1 - obj.e_LSO.^2)/(6 + 2.*obj.e_LSO)).^(3/2) ;
                    % x0: Vector with initial conditions for the ODE system. Here is defined to be equal to the LSO.
                    x0 = [obj.Phi_LSO, nu_LSO, obj.e_LSO, obj.Gamma_LSO, obj.alpha_LSO ] ; 
                    t_min   = 0 ;                              % final integration time.
                    dt      = 1/(10.*nu_LSO);
                    points  = (floor(1/dt.*abs(obj.t_max-t_min)));  % grid points for integration. Defines resolution of waveform
                    %t_span  = linspace(t_min,  obj.t_max, points) ;     % time limits for integration
                    t_span  = linspace(obj.t_max,  t_min, points);      % time limits for integration
                    rtol    = 10.^-10 ;                                % relative tolerance for ODE solver
                    atol    = 10.^-10   ;                              % absolute tolerance for ODE solver
                    
                    [T1,sol] = ode45(@(t,y) eqs(t,y,obj),t_span,x0);
                   
                    t = flip(T1); 
                    Phi=sol(:,1); Phi = flip(Phi); 
                    nu=sol(:,2);  nu  = flip(nu);
                    e=sol(:,3);  e   = flip(e); 
                    gamma=sol(:,4); gamma = flip(gamma); 
                    alpha=sol(:,5); alpha = flip(alpha); 
           end
           
           function [dy]=eqs(t,y,obj)
                % eqs()
                % Note: The initial values for the variables (Phi, y(5), gamma, e) at LSO are functionined in global_parameters.py.
                %     """
                %     functioninition of system of kludge orbital equations as in Barrack&Cutler2004.
                %     Extend the function arguments in case needed to include other variables (eccentricity, angles etc.)
                % 
                % Parameters
                % ----------
                % t : 1D array
                % y : 1D array
                % M0: central mass M0 (solar mass)
                % u0: orbiting mass M0 (solar mass)
                % spin0: spin of central black hole (G./c)
                % 
                % Returns
                % -------
                % the l.h.s of the system of ODE equations
                %      """
                % Vector of orbital variables: Phase, frequency, angles gamm and alph
                %     fai= y(1);
                %     nu = y(2);
                %     e = y(3);
                %     gamma = y(4) ;
                %     alpha = y(5) ;

                %x    = ((2*Pi*Parameters.G*M*y(2))./Parameters.c.^3);  % Post-Newtonian parameter
                M    = obj.M0*Common_CF.Ms  ;    % Restore Parameters in central mass
                mu    = obj.mu0*Common_CF.Ms   ;   % Restore Parameters in orbiting mass 
                spin = obj.spin0*(Common_CF.G/Common_CF.c);% Restore Parameters in spin

                %%%%%% System of ODEs for {?, nu, e, ?, ?} %%%%%%
                dy = zeros(5,1);
                dy(1) = 2*pi*y(2);
                dy(2) = ((96./(10*pi))*((Common_CF.c.^6./Common_CF.G.^2)*mu./M.^3)*(((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3).^(11./3))*(1 - y(3).^2).^(-9./2)*((1 - y(3).^2)*(1 + (73./24)*y(3).^2 + (37./96)*y(3).^4) + ((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3).^(2./3)*((1273./336) - (2561./224)*y(3).^2 - (3885./128)*y(3).^4 - (13147./5376)*y(3).^6) - ((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3)*(Common_CF.c./Common_CF.G*spin)*cos(obj.lambda)*(1 - y(3).^2).^(-1./2)*((73./12) + (1211./24)*y(3).^2 + (3143./96)*y(3).^4 + (65./64)*y(3).^6))); 
                dy(3) = (-(y(3)./15)*((Common_CF.c.^3./Common_CF.G)*mu./M.^2)*(((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3).^(8./3))*((1 - y(3).^2).^(-7./2))*( (1 - y(3).^2)*(304 + 121*y(3).^2)*(1+12*((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3).^(2/3)) - 1./56*((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3).^(2/3)*( 8*16705 + 12*9082*y(3).^2 - 25211*y(3).^4)) + y(3)*((Common_CF.c.^3./Common_CF.G)*mu./M.^2)*(Common_CF.c./Common_CF.G*spin)*cos(obj.lambda)* ((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3).^(11./3)*(1-y(3).^2).^(-4)*( ( 1364/5) + (5032/15)*y(3).^2 + (263/10)*y(3).^4));
                dy(4) = 6*pi*y(2)*((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3).^(2./3)*((1 - y(3).^2).^(-1))* (1 - 2*(((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3).^(1./3))*(Common_CF.c./Common_CF.G*spin)*cos(obj.lambda)*(1 - y(3).^2).^(-1./2)+ (1./4)*((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3).^(2./3)*(1 - y(3).^2).^(-1)*(26 - 15*y(3).^2));
                dy(5) = 4*pi*y(2)*((2*pi*Common_CF.G*M*y(2))./Common_CF.c.^3)*(Common_CF.c./Common_CF.G*spin)*(1 - y(3).^2).^(-3./2);
           end
             
           function [Aplus_n,Across_n]=get_Aplus_A_cross(obj,Ln, a_n,b_n, c_n, ytilde) 
                   Aplus_n = -(1.0 + Ln.^2) .* ( a_n .* cos(2.0 .* ytilde) - b_n .* sin(2.0 .* ytilde) ) + (1.0 - Ln.^2) .* c_n;
                   Across_n = 2.0 .* Ln .* (b_n .* cos(2.0 .* ytilde) + a_n .* sin(2.0 .* ytilde));
           end 
           
            function [a_n,b_n,c_n]=get_an_bn_cn(obj,A,n,e,Fi)
                   a_n = -n .* A .* ( besselj(n-2,n*e) - 2 .* e .* besselj(n-1,n*e) + (2.0/n) .* besselj(n,n*e) + 2 .* e .* besselj(n+1,n*e) - besselj(n+2,n*e) ) .* cos(n.*Fi);
                   b_n = -n .* A .* ((1 - e.^2).^0.5) .* ( besselj(n-2,n*e) - 2 .* besselj(n,n*e) + besselj(n+2,n*e) ) .* sin(n.*Fi);
                   c_n = 2 .* A .* besselj(n,n*e) .* cos(n.*Fi);
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
                    ret=Common_CF.c./Common_CF.H0.*(1+z).*sum(1./sqrt(O_M.*(1+zz).^3+O_A)).*dzz; %Bennet C L et al 2003 Astrophys. J. S. 148 1
        end
       
          function [fsg,hsg]=Fourier_tran(obj,t,h)
                dt = abs(t(2) - t(1));
                Dt=max(t)-min(t);
                pp=fft(h)*dt;
                w1=pi./dt;
                dw=2*pi/Dt;
                w=-w1:dw:w1;
                fff=w./(2.*pi)./(1+ obj.Zf);
                sdf=(fftshift(pp'));
                N = ceil(length(fff)/2)+1;
                fsg=fff(N:length(fff));
                hsg = sdf(N:length(fff))';
          end
       end
   end
   