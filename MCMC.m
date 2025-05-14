classdef MCMC
    % Define the MCMC class with properties and methods. 
    properties (Constant = true)  
        Ntotal         % Total number of iterations for the MCMC run.
        printerval     % Interval at which to print status updates.
    end
    methods

       function [samples,lps]=MCMC_run(obj,data_f,freq_bin,true_vals,cov_Matrix,variances,params_low,params_sh,params_high,Ntotal,printerval)
            % MCMC_RUN - Runs the MCMC algorithm.
            %   obj: instance of the MCMC class.
            %   data_f: observed data in the frequency domain.
            %   freq_bin: frequency bins.
            %   true_vals: true parameter values (initial guess).
            %   cov_Matrix: covariance matrix.
            %   variances: noise variances.
            %   params_low: lower bounds for the parameters.
            %   params_high: upper bounds for the parameters.
            %   Ntotal: total number of iterations.
            %   printerval: interval for printing status.
            
            %rand('seed',2); % Set the seed for reproducibility.
            % Initialize the parameters array with the first row being the true values.
            paramst = [true_vals params_sh];
            % Generate the initial waveform in the frequency domain.
            WF = Waveform(paramst); % Assuming 'waveform' is a class that generates waveforms.
            signal_init_f = WF.waveform_fd(freq_bin); % Generate the initial waveform.
            % Compute the periodogram using the initial signal.
            pdgrm = abs(data_f - signal_init_f).^2;

            cov(1,:,:)= cov_Matrix;
            dim = size(true_vals,2);
            chains = 2*dim+2;
            % start parpool by codes.
            delete(gcp('nocreate'));         % stop the process before start a new run.
            numCore = chains+1;   % get the maxmium core num of PC.
            parpool(numCore-1);              % start parpool.

            %chains = 1;
            accept_reject_count=0;
            parfor jj=1:chains
                %rand('seed',2*jj); % Set the seed for reproducibility.
                lp     = 0;
                sample = true_vals; % Initialize the log posterior array.
                lp_store = [lpost(obj,pdgrm,variances,true_vals,params_low,params_high)];  
                cov_Tem = reshape(cov, dim, dim);
            for ii =2:Ntotal  
                if mod(ii,printerval) == 0
                    str = ['ii= ', num2str(ii)];
                    disp(str); 
                end
                % Retrieve the previous log posterior and parameter values.
                lp_prev = lp_store;
                prev_vec = sample(ii-1,:);

                
                params_prop= mvnrnd(prev_vec, (1/4).*cov_Tem);%the mean valua is prev_vec,and the standard covariance is Cov_Matrix
                if (mod(ii,50)==0)
                    cov_Tem = updateCovarianceMatrix(obj,sample,cov_Tem);
                end                                            
                % Generate the proposed signal.
                % params_prop= mvnrnd(prev_vec, (1/4).*cov_Matrix);%
                paramst = [params_prop params_sh];
                WF = Waveform(paramst);
                signal_prop_f = WF.waveform_fd(freq_bin);
                pdgrm_prop = abs(data_f - signal_prop_f).^2;
                
                % Compute the log posterior for the proposed parameters.
                lp_prop = lpost(obj,pdgrm_prop,variances,params_prop,params_low,params_high);
                % Decide whether to accept or reject the proposal based on the Metropolis-Hastings criterion.
               
                 if accept_reject(obj,lp_prop,lp_prev) == 1
                    sample(ii,:) = params_prop;
                    accept_reject_count = [accept_reject_count;1];
                    lp_store = lp_prop;
                else
                    % Reject the proposal.
                    sample(ii,:) = sample(ii-1,:);
                    accept_reject_count = [accept_reject_count;0];
                 end
                lp = [lp;lp_store];
            end
            lp_storej(jj)=lp_store;
            samples(jj,:,:) = sample(:,:);
            lps(jj,:)=lp';
            end
       end

      function [sam_con,lp_con]=MCMC_contin(obj,data_f,freq_bin,cov_Matrix,variances,params_low,params_sh,params_high,samples,lps,Ntotal,accept_reject_count,printerval)
          
            delta_f = freq_bin(2) - freq_bin(1);
            Nini = size(samples,2);
            dim = size(samples,3);
            chains = size(samples,1);
            lpss = [];
            sampless=[];
     
            % start parpool by codes.
            delete(gcp('nocreate'));         % stop the process before start a new run.
            numCore = chains;   % get the maxmium core num of PC.
            parpool(numCore-1);              % start parpool.
            cov(1,:,:)= cov_Matrix;
            parfor jj = 1:chains
                   temp_lp = reshape(lps(jj,:),Nini,1);
                   temp_chains = reshape(samples(jj,:,:),Nini,dim);
                   cov_Tem = reshape(cov, dim, dim);
                   
            for ii = Nini+1:Ntotal
                % Print status updates at specified intervals.
                if mod(ii,printerval) == 0
                    str = ['ii= ', num2str(ii)];
                    disp(str);
                end
                % Retrieve the previous log posterior and parameter values.
                prev_vec = temp_chains(ii-1,:);
                % Propose new parameter values from a multivariate normal distribution.
                %params_prop = mvnrnd(prev_vec, (1/4)*cov_Matrix);
             
                params_prop= mvnrnd(prev_vec, (1/4).*cov_Tem);%the mean valua is prev_vec,and the standard covariance is Cov_Matrix
                if (mod(ii,50)==0)
                    cov_Tem = updateCovarianceMatrix(obj,temp_chains,cov_Tem);
                end
                lp_prev = temp_lp(ii-1);
                % Generate the proposed signal.
                paramst = [params_prop params_sh];
                WF = Waveform(paramst);
                signal_prop_f = WF.waveform_fd(freq_bin);
                %signal_prop_f = WF.waveform_td(sta_t);
                % Compute the periodogram for the proposed signal.
                pdgrm_prop = abs(data_f - signal_prop_f).^2;
                
                % Compute the log posterior for the proposed parameters.
                lp_prop = lpost(obj,pdgrm_prop,variances,params_prop,params_low,params_high);
                % Decide whether to accept or reject the proposal based on the Metropolis-Hastings criterion.
                if accept_reject(obj,lp_prop,lp_prev) == 1
                %if MCMC.accept_reject(lp_prop,lp_prev) == 1
                    % Accept the proposal.
                    temp_chains(ii,:)=params_prop;
                    %accept_reject_count = [accept_reject_count;1];
                    temp_lp(ii)=lp_prop;
                else
                    % Reject the proposal.
                    temp_chains(ii,:)=prev_vec;
                    %accept_reject_count = [accept_reject_count;0];
                    temp_lp(ii)=lp_prev;
                end
            end
                sampless = [sampless;temp_chains];
                lpss     = [lpss;temp_lp]
            end
            sam_con = shape_down(obj,sampless,chains,Ntotal);
            lp_con = shape_down(obj,lpss,chains,Ntotal);
           end
      
   function [samples_tol,lps_tol,R] = converg(obj, sam_con, lp_con)
            Ntotal = size(sam_con,2);
            n = floor(Ntotal./2./2);
            chains = size(sam_con,1);
               samplesMeff1(:,:,:) = sam_con(:,(2*n+1):3*n,:);
               samplesMeff2(:,:,:) = sam_con(:,(3*n+1):4*n,:);
               samplesMeff=[ samplesMeff1; samplesMeff2];
               lpMeff1(:,:) = lp_con(:,(2*n+1):3*n);
               lpMeff2(:,:) = lp_con(:,(3*n+1):4*n);
               lpMeff=[lpMeff1;lpMeff2];
               m=2*chains;
               for jj=1:m
                   samples_u(jj,:)=sum(samplesMeff(jj,:,:))./n;   
               end
               mean = sum(samples_u(:,:))./m;

               B = sum((samples_u(:,:)-mean).^2./(m-1));

               W = 0;
               for jj=1:m  
               samples_j(:,:)=samplesMeff(jj,:,:);
               s2 = sum((samples_j(:,:) - ones(n,size(sam_con,3)).*samples_u(jj,:) ).^2)./(n-1);
               W  = (W + s2);
               end
               W  = W./m;
               sigma = (n-1)./n.*W+B;
               R = sqrt(sigma./W);
               if any(R > 1.05)
                   disp('sampling points may be insufficient£¡');
               end
               samples_tol = reshape(samplesMeff,[m*n,size(sam_con,3)]); 
               lps_tol = reshape(lpMeff,[m*n,1]); 
               end
        
        function [res]=lpost(obj,pdgrm,variances,params,params_low,params_high)
            % LPOST - Computes the log posterior probability.
            %   pdgrm: periodogram.
            %   variances: noise variances.
            %   params: parameter vector.
            %   params_low: lower bounds for the parameters.
            %   params_high: upper bounds for the parameters.
            
            % Compute the prior probabilities for each parameter.
            re = zeros(size(params));
            for ii=1:length(params)
                re(ii) = lprior(obj,params(ii),params_low(ii),params_high(ii));
            end
            
            % Compute the log likelihood.
            res = sum(re) + llike(obj,pdgrm,variances);
        end
        
        function [re]=lprior(obj,M_chirp,M_chirp_low,M_chirp_high)
            % LPRIOR - Computes the log prior probability for a single parameter.
            %   M_chirp: parameter value.
            %   M_chirp_low: lower bound for the parameter.
            %   M_chirp_high: upper bound for the parameter.  
            % Check if the parameter is within the valid range.
            if (M_chirp < M_chirp_low) || (M_chirp > M_chirp_high)
                re = -1e100; % Assign a very low probability if out of bounds.
                % str = ['rejected '];
                % disp(str);
            else
                re = 0; % Uniform prior.
            end
        end
        
        function [re]=llike(obj,pdgrm,variances)
            % LLIKE - Computes the log likelihood.
            %   pdgrm: periodogram.
            %   variances: noise variances.
            
            % Compute the log likelihood using the Whittle likelihood.
            re = -0.5 * sum(pdgrm ./ variances);
        end
        
        function [re]=accept_rate(obj,parameter)
            % ACCEPT_RATE - Computes the acceptance rate for a sequence of parameter samples.
            %   parameter: sequence of samples of a parameter.
            
            % Initialize the rejection counter.
            rejections = 0;
            for i = 1:(length(parameter)-1)
                % Count the number of rejections.
                rejections = rejections + (parameter(i + 1) == parameter(i));
            end
            reject_rate = rejections / (length(parameter) - 1); % Compute the rejection rate.
            re = 1 - reject_rate; % Compute the acceptance rate.
        end
        
        function [re]=accept_reject(obj,lp_prop,lp_prev)
            % ACCEPT_REJECT - Computes the decision to accept or reject a proposal.
            %   lp_prop: proposed log posterior.
            %   lp_prev: previous log posterior.
            
            % Generate a random number between 0 and 1.
            u = unifrnd(0,1);
            
            % Compute the log acceptance probability.
            r = min(0, lp_prop-lp_prev);
            
            % Decide whether to accept (1) or reject (0) the proposal.
            if log(u) < r
                re = 1; % Accept the proposal.
            else
                re = 0; % Reject the proposal.
            end
        end

        
        function [sa]=shape_down(obj,sampless,chains,Ntotal)
                for ii=1:chains
                    sa(ii,:,:) = sampless((ii-1)*Ntotal+1:ii*Ntotal,:);
                end
        end

        function Cov_Matrix = updateCovarianceMatrix(obj,samples,Cov_Matrix)
            % Number of parameters
            N = size(samples, 2); % Number of parameters
            lambda = size(samples, 1); % Number of samples
            mu = floor(lambda / 2); % Number of best samples
            weights = log(mu + 0.5) - log(1:mu); % Weights for the best samples
            weights = weights / sum(weights); % Normalize weights
            mueff = 1 / sum(weights .^ 2); % Variance-effective size of mu
            % Compute the mean of the samples
            xmean = mean(samples, 1);
             ps = zeros(N, 1);
            pc = zeros(N, 1);
            invsqrtC = eye(N); % Inverse square root of the covariance matrix
            sigma = 1.0; % Step size
            cs = 1 / sqrt(N + 1); % Time constant for cumulation for step-size control
            cc = 2 / (N + 2); % Time constant for cumulation for the covariance matrix
            c1 = 2 / ((N + 1.3)^2 + mueff); % Learning rate for rank-one update
            cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((N + 2)^2 + mueff)); % Learning rate for rank-mu update
            % Compute the evolution paths
            xold = samples(1, :);
            ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * invsqrtC * (xmean - xold)' / sigma;
            hsig = sum(ps.^2) / (1 - (1 - cs)^(2 * lambda / mu)) / N < 2 + 4 / (N + 1);
            pc = (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mueff) * (xmean - xold)' / sigma;     
            % Adapt covariance matrix C
            artmp = (1 / sigma) * (samples(1:mu, :) - repmat(xold, mu, 1)); % mu difference vectors
            C = (1 - c1 - cmu) * Cov_Matrix + c1 * (pc * pc' + (1 - hsig) * cc * (2 - cc) * Cov_Matrix) + cmu * artmp' * diag(weights) * artmp;
            % Update the inverse square root of the covariance matrix
            [V, D] = eig(C);
            invsqrtC = V * diag(1 ./ sqrt(diag(D))) * V';
            % Return the updated covariance matrix
            Cov_Matrix = C;
        end
    end
end
