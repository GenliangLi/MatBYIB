classdef Fisher_Matrix
    %FISHER_MATRIX  Class for computing Fisher information matrix and parameter covariance
    %   This class implements methods for:
    %       - Fisher matrix calculation via numerical derivatives
    %       - Parameter covariance estimation
    %       - Frequency-domain inner products for gravitational wave analysis
    %
    %   Example:
    %       fisher = Fisher_Matrix();
    %       covMat = fisher.Matrix(params, f, PSD, nParams);
    
    methods
        function [cov_Matrix] = Matrix(obj, params, freq_bin, PSD, NN)
            %MATRIX  Compute parameter covariance matrix from Fisher information matrix
            %
            %   [cov_Matrix] = Matrix(obj, params, freq_bin, PSD, NN)
            %
            %   Inputs:
            %       params    - Parameter vector [1×N] (physical units)
            %       freq_bin  - Frequency bins [1×M] (Hz)
            %       PSD       - Power spectral density [1×M] (strain^2/Hz)
            %       NN        - Number of parameters (scalar)
            %
            %   Outputs:
            %       cov_Matrix - Parameter covariance matrix [N×N]
            %
            %   Algorithm:
            %       1. Computes partial derivatives via central differences
            %       2. Constructs Fisher matrix using noise-weighted inner products
            %       3. Inverts Fisher matrix to obtain covariance
            
            % Initialize derivative matrix [Nparams × Nfreq]
            diff_vec = zeros(length(params), length(freq_bin));
            
            % Compute numerical derivatives for each parameter
            for ii = 1:NN
                diff_vec(ii,:) = diff_param(obj, ii, params, freq_bin);
            end
            diff_vec(isnan(diff_vec)) = 0;

            % Initialize Fisher matrix
            fish_mix = eye(NN);  % Preallocate for speed
            delta_f = freq_bin(2) - freq_bin(1);  % Frequency resolution
            
            % Build Fisher matrix
            for ii = 1:NN
                for jj = 1:NN
                    fish_mix(ii,jj) = inner_prod(obj, diff_vec(ii,:), ...
                                             diff_vec(jj,:), PSD, delta_f);
                end
            end
            % Compute covariance matrix (matrix inverse of Fisher matrix)
            cov_Matrix = inv(fish_mix);
        end
        
        function [diff_vec] = diff_param(obj, ii, params, freq_bin)
            %DIFF_PARAM  Compute numerical derivative via central differences
            %
            %   [diff_vec] = diff_param(obj, ii, params, freq_bin)
            %
            %   Inputs:
            %       ii       - Parameter index to differentiate
            %       params   - Full parameter vector
            %       freq_bin - Frequency bins [1×M] (Hz)
            %
            %   Outputs:
            %       diff_vec - Frequency-domain derivative [1×M]
            %
            %   Notes:
            %       - Uses 2nd-order accurate central difference scheme
            %       - Step size scales with parameter magnitude
            
            % Set adaptive step size (1e-10 of parameter value)
            delta = 1e-10 * params(ii);
            
            % Create parameter perturbations
            params_p = params;
            params_m = params;
            params_p(ii) = params(ii) + 0.5*delta;
            params_m(ii) = params(ii) - 0.5*delta;
            
            % Compute perturbed waveforms
            wave_p = waveform(params_p);
            hf_p = wave_p.waveform_fd(freq_bin);
            
            wave_m = waveform(params_m);
            hf_m = wave_m.waveform_fd(freq_bin);
            
            % Central difference derivative
            diff_vec = (hf_p - hf_m) / delta;
        end
        
        function [re] = inner_prod(obj, sig1_f, sig2_f, PSD, delta_f)
            %INNER_PROD  Noise-weighted inner product in frequency domain
            %
            %   [re] = inner_prod(obj, sig1_f, sig2_f, PSD, delta_f)
            %
            %   Inputs:
            %       sig1_f  - First signal [1×M] (complex)
            %       sig2_f  - Second signal [1×M] (complex)
            %       PSD     - Power spectral density [1×M] (strain^2/Hz)
            %       delta_f - Frequency bin width (Hz)
            %
            %   Outputs:
            %       re      - Inner product result (real scalar)
            %
            %   Formula:
            %       re = 4Δf·Re[∑(sig1(f)·sig2*(f)/PSD(f))]
            
            re = (4 * delta_f) * real(sum(sig1_f .* conj(sig2_f) ./ PSD));
        end
    end
end
