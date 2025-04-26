classdef plotFigure
    %PLOTFIGURE Visualization class for gravitational wave analysis results
    %
    %   This class provides methods for plotting:
    %   - Time-domain and frequency-domain waveforms
    %   - Corner plots for parameter distributions
    %   - Posterior distributions
    %   - Orbital evolution plots
    %
    %   Methods:
    %     plot_t_fig     - Plot time-domain waveforms
    %     plot_f_fig     - Plot frequency-domain waveforms
    %     cornerplot     - Create corner plots for parameter distributions
    %     Plot_Posterior - Plot posterior distributions
    %     plotOrbit      - Plot orbital evolution parameters
    properties 
    end
    methods
        function plot_t_fig(obj,t,ht,ty)
            %PLOT_T_FIG Plot time-domain waveforms
            %
            % Inputs:
            %   t  : Time vector [Nx1] (seconds)
            %   ht : Waveform amplitudes [Nx1]
            %   ty : Plot type (1=loglog, 2=linear)
            figure
            if ty==1
                loglog(t,ht,'k');hold on
                set(gca,'fontsize',12);
                set(gca,'LineWidth',1); 
                set(get(gca,'xlabel'),'Fontsize',14);
                set(get(gca,'ylabel'),'Fontsize',14);
                xlabel('Time(s)','FontSize',14,'Fontname','ËÎÌå'); 
                ylabel('h_t','FontSize',14,'Fontname','ËÎÌå'); 
            else
                plot(t,ht,'k');hold on
                set(gca,'fontsize',12);
                set(gca,'LineWidth',1); 
                set(get(gca,'xlabel'),'Fontsize',14);
                set(get(gca,'ylabel'),'Fontsize',14);
                xlabel('Time(s)','FontSize',14,'Fontname','ËÎÌå'); 
                ylabel('h_t','FontSize',14,'Fontname','ËÎÌå'); 
            end
        end
        function plot_f_fig(obj,f,hf,ty)
            %PLOT_F_FIG Plot frequency-domain waveforms
            %
            % Inputs:
            %   f  : Frequency vector [Nx1] (Hz)
            %   hf : Spectral amplitudes [Nx1]
            %   ty : Plot type (1=loglog, 2=linear)
            figure
            if ty==1
                loglog(f,hf,'k');hold on
                set(gca,'fontsize',12);
                set(gca,'LineWidth',1); 
                set(get(gca,'xlabel'),'Fontsize',14);
                set(get(gca,'ylabel'),'Fontsize',14);
                xlabel('Frequency(Hz)','FontSize',14,'Fontname','ËÎÌå'); 
                ylabel('h_f','FontSize',14,'Fontname','ËÎÌå'); 
            else
                plot(f,hf,'k');hold on
                set(gca,'fontsize',12);
                set(gca,'LineWidth',1); 
                set(get(gca,'xlabel'),'Fontsize',14);
                set(get(gca,'ylabel'),'Fontsize',14);
                xlabel('Frequency(Hz)','FontSize',14,'Fontname','ËÎÌå'); 
                ylabel('h_f','FontSize',14,'Fontname','ËÎÌå'); 
            end
        end
        function [fig,ax]=cornerplot(obj,data, varargin)
            %CORNERPLOT Create corner plot for parameter distributions
            %
            % Inputs:
            %   data    : Parameter samples [Nsamples x Nparams]
            %   varargin: Optional name-value pairs
            %     - names  : Cell array of parameter names
            %     - truths : True parameter values for reference
            %     - titles : Plot titles
            %
            % Outputs:
            %   fig : Figure handle
            %   ax  : Axes handles
                    format long            
                    if length(size(data)) ~= 2
                        error('x must be 2D.')
                    end       
                    nDims = min(size(data));  
                    if nDims ~= size(data,2)
                        data = data';
                    end      
                    names = {};
                    truths = [];
                    bounds = [];
                    bounds_supplied = true;
                    top_margin = 0;        
                    gutter = [.004 .004];
                    margins = [.1 .01 .12 .1];         
                    if nargin > 1
                        names = varargin{1}
                        if ~isempty(names) && ~(iscell(names) && length(names) == nDims)
                            error('NAMES must be a cell array with length equal to the number of dimensions in your data.')
                        end
                        if nargin > 2
                            truths = varargin{2}
                            if ~isempty(truths) && ~(isfloat(truths) && numel(truths) == nDims)
                                error('TRUTHS must be a vector with length equal to the number of dimensions in your data.')
                            end
                            if nargin > 3
                                    titles=varargin{3};
                            end
                            end
                    end    
                    if isempty(bounds) | all(bounds==0)
                        bounds = nan(2,nDims);
                        bounds_supplied = false;
                    end
                    fig = figure;
                    set(gcf, 'color', 'w')
                    ax = nan(nDims+top_margin,nDims);
                    hist_bins = 30;
                    lines = 10;
                    res = 2^5;
                    linewidth = 1;
                    axes_defaults = struct('tickdirmode','manual',...
                        'tickdir','in',...
                        'ticklength',[.035 .035],...
                        'box','on',...
                        'xticklabel',[],...
                        'yticklabel',[]); 
                    std_devs = std(data);   
                    for i = 1:nDims
                        if ~bounds_supplied
                            bounds(:,i) = [min(data(:,i)) max(data(:,i))];
                        end
                        
                        for t = 1:top_margin
                            ax(i-1+t,i) = tight_subplot(obj,nDims+top_margin, nDims, i-1+t, i, gutter, margins);
                            set(gca,'visible','off','xlim',bounds(:,i));
                        end
                    
                        truncated_data = data;
                        truncated_data(truncated_data(:,i)<bounds(1,i) | truncated_data(:,i)>bounds(2,i),i) = nan;
                        
                        ax(i+top_margin,i) = tight_subplot(obj,nDims+top_margin, nDims, i+top_margin, i, gutter, margins);
                        edges=linspace(min(truncated_data(:,i)),max(truncated_data(:,i)),hist_bins);
                        h=histogram(truncated_data(:,i), 'BinEdges', edges,'normalization', 'probability', 'displaystyle', 'stairs', 'edgecolor', 'k')
                        set(gca,'xlim', bounds(:,i),'ylim', [0 max(h.Values)], axes_defaults,'ytick',[]);
                        
                        x_ticks = get(gca,'XTick');
                        x_tick_labels = arrayfun(@num2sci_str, x_ticks, 'UniformOutput', false);
                        set(gca,'XTickLabel', x_tick_labels);
                    
                        hold on;
                        mean_val = mean(data(:,i));
                        plot([mean_val+std_devs(i) mean_val+std_devs(i)], [0 max(h.Values)], 'r--', 'linewidth',linewidth);hold on
                        plot([mean_val-std_devs(i) mean_val-std_devs(i)], [0 max(h.Values)], 'r--', 'linewidth',linewidth);hold on
                        title(titles{i})
                        if i == nDims
                            set(gca,'xticklabelmode','auto')
                        end
                        
                        if ~isempty(truths)
                            hold on
                            plot([truths(i) truths(i)], [0 1], 'r--', 'linewidth',linewidth)
                        end
                      
                        if ~isempty(names)
                            if i == 1
                                ylabel(names{i});
                            end
                            if i == nDims
                                xlabel(names{i})
                            end
                        end
                        
                    end
                    
                    if nDims > 1
                        for d1 = 1:nDims-1
                            for d2 = d1+1:nDims
                                [~, density, X, Y] = kde2d(obj,[data(:,d1) data(:,d2)],res,[bounds(1,d1) bounds(1,d2)],[bounds(2,d1) bounds(2,d2)]);
                    
                                ax(d2+top_margin,d1) = tight_subplot(obj,nDims+top_margin, nDims, d2+top_margin, d1, gutter, margins);
                                contour(X,Y,10^10*density, lines)
                                
                                set(gca,'xlim',bounds(:,d1),'ylim',bounds(:,d2), axes_defaults);
                                
                                if ~isempty(truths)
                                    yl = get(gca,'ylim');
                                    xl = get(gca,'xlim');
                                    hold on
                                    plot(xl, [truths(d2) truths(d2)],'r--', 'linewidth',linewidth)
                                    plot([truths(d1) truths(d1)], yl,'r--', 'linewidth',linewidth)
                                end
                                x_ticks = get(gca,'XTick');
                                x_tick_labels = arrayfun(@num2sci_str, x_ticks, 'UniformOutput', false);
                                set(gca,'XTickLabel', x_tick_labels);
                        
                                if d1 == 1
                                    if ~isempty(names)
                                        ylabel(names{d2})
                                    end
                                    set(gca,'yticklabelmode','auto')
                                end
                                if d2 == nDims
                                    if ~isempty(names)
                                        xlabel(names{d1})
                                    end
                                    set(gca,'xticklabelmode','auto')
                                end
                            end
                            
                            row = ax(1+top_margin+d1,:);
                            row = row(~isnan(row));
                            row = row(1:d1);
                            
                            col = ax(:,d1);
                            col = col(~isnan(col));
                            col = col(1:end);
                            
                            linkaxes(row, 'y');
                            linkaxes(col, 'x');
                        end
                    end
            end

            function [bandwidth,density,X,Y]=kde2d(obj,data,n,MIN_XY,MAX_XY)
                %KDE2D 2D Kernel Density Estimation
                %
                % Inputs/Outputs match original implementation
                    global N A2 I
                    if nargin<2
                        n=2^8;
                    end
                    n=2^ceil(log2(n));
                    N=size(data,1);
                    if nargin<3
                        MAX=max(data,[],1); MIN=min(data,[],1); Range=MAX-MIN;
                        MAX_XY=MAX+Range/4; MIN_XY=MIN-Range/4;
                    end
                    scaling=MAX_XY-MIN_XY;
                    if N<=size(data,2)
                        error('data has to be an N by 2 array where each row represents a two dimensional observation')
                    end
                    transformed_data=(data-repmat(MIN_XY,N,1))./repmat(scaling,N,1);
                    initial_data=ndhist(obj,transformed_data,n);
                    a= dct2d(obj,(initial_data));
                      I=(0:n-1).^2; A2=a.^2;

                     t_star=fzero( @(t)(t-evolve(obj,t)),[0,0.1]);

                    p_02=func(obj,[0,2],t_star);p_20=func(obj,[2,0],t_star); p_11=func(obj,[1,1],t_star);
                    t_y=(p_02^(3/4)/(4*pi*N*p_20^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
                    t_x=(p_20^(3/4)/(4*pi*N*p_02^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
                    a_t=exp(-(0:n-1)'.^2*pi^2*t_x/2)*exp(-(0:n-1).^2*pi^2*t_y/2).*a; 
                    if nargout>1
                        density=idct2d(obj,a_t)*(numel(a_t)/prod(scaling));
                        [X,Y]=meshgrid(MIN_XY(1):scaling(1)/(n-1):MAX_XY(1),MIN_XY(2):scaling(2)/(n-1):MAX_XY(2));
                    end
                    bandwidth=sqrt([t_x,t_y]).*scaling; 
            end
    
            function binned_data=ndhist(obj,data,M)
                    [nrows,ncols]=size(data);
                    bins=zeros(nrows,ncols);
                    for i=1:ncols
                        [dum,bins(:,i)] = histc(data(:,i),[0:1/M:1],1);
                        bins(:,i) = min(bins(:,i),M);
                    end
                    binned_data = accumarray(bins(all(bins>0,2),:),1/nrows,M(ones(1,ncols)));
            end

            function out=K(obj,s)
                     out=(-1)^s*prod((1:2:2*s-1))/sqrt(2*pi);
            end
        
            function  [out,time]=evolve(obj,t)
                    global N
                    Sum_func = func(obj,[0,2],t) + func(obj,[2,0],t) + 2*func(obj,[1,1],t);
                    time=(2*pi*N*Sum_func)^(-1/3);
                    out=(t-time)/time;
            end

            function out=func(obj,s,t)
                        global N
                        if sum(s)<=4
                            Sum_func=func(obj,[s(1)+1,s(2)],t)+func(obj,[s(1),s(2)+1],t); const=(1+1/2^(sum(s)+1))/3;
                            time=(-2*const*K(obj,s(1))*K(obj,s(2))/N/Sum_func)^(1/(2+sum(s)));
                            out=psis(obj,s,time);
                        else
                            out=psis(obj,s,t);
                        end
            end
                     
            function out=psis(obj,s,Time)
                        global I A2
                        w=exp(-I*pi^2*Time).*[1,.5*ones(1,length(I)-1)];
                        wx=w.*(I.^s(1));
                        wy=w.*(I.^s(2));
                        out=(-1)^sum(s)*(wy*A2*wx')*pi^(2*sum(s));
            end

            function data=dct2d(obj,data)
                    [nrows,ncols]= size(data);
                    if nrows~=ncols
                        error('data is not a square array!')
                    end
                    w = [1;2*(exp(-i*(1:nrows-1)*pi/(2*nrows))).'];
                    weight=w(:,ones(1,ncols));
                    data=dct1d(dct1d(data)')';
                        function transform1d=dct1d(x)
                            x = [ x(1:2:end,:); x(end:-2:2,:) ];
                            transform1d = real(weight.* fft(x));
                        end
            end

            function data = idct2d(obj,data)
                    [nrows,ncols]=size(data);
                    w = exp(i*(0:nrows-1)*pi/(2*nrows)).';
                    weights=w(:,ones(1,ncols));
                    data=idct1d(idct1d(data)');
                        function out=idct1d(x)
                            y = real(ifft(weights.*x));
                            out = zeros(nrows,ncols);
                            out(1:2:nrows,:) = y(1:nrows/2,:);
                            out(2:2:nrows,:) = y(nrows:-1:nrows/2+1,:);
                        end
            end
            function Plot_Posterior(obj)
                    figure
                    subplot(2,3,1)
                    hist(samples(:,1), 30, 'FaceColor', 'b'); hold on
                    line([(M0) M0], [0 max(histcounts(samples(:,1), 30))./2], 'Color', 'r'); hold on
                    ylabel('P({M}_{0}|m_0,z,f_0,e_0,\phi_0)')
                    xlabel('{M}_{0}')
                    title('Center Mass - {M}_{0}')
                    legend('Posterior', 'True Value')

                    subplot(2,3,2)
                    hist(samples(:,2), 30, 'FaceColor', 'b'); hold on
                    line([m0 m0], [0 max(histcounts(samples(:,2), 30))./2], 'Color', 'r'); hold on
                    ylabel('P(\mu|M_0,z,f_0,e_0,\phi_0)')
                    xlabel('\mu')
                    title('Small Mass - \mu')
                    legend('Posterior', 'True Value')

                    subplot(2,3,3)
                    hist(samples(:,3), 30, 'FaceColor', 'b'); hold on
                    line([Zf Zf], [0 max(histcounts(samples(:,3), 30))./2], 'Color', 'r'); hold on
                    ylabel('P(z|M_0,m_0,f_0,e_0,\phi_0)')
                    xlabel('z')
                    title('Redshift - z')
                    legend('Posterior', 'True Value')

                    subplot(2,3,4)
                    hist(samples(:,4), 30, 'FaceColor', 'b'); hold on
                    line([forb0 forb0], [0 max(histcounts(samples(:,4), 30))./2], 'Color', 'r'); hold on
                    ylabel('P(f0|M_0,m_0,z,e_0,\phi_0)')
                    xlabel('f0')
                    title('Frequency - f0')
                    legend('Posterior', 'True Value')

                    subplot(2,3,5)
                    hist(samples(:,5), 30, 'FaceColor', 'b'); hold on
                    line([e0 e0], [0 max(histcounts(samples(:,5), 30))./2], 'Color', 'r'); hold on
                    ylabel('P(e0|M_0,m_0,z,f_0,\phi_0)')
                    xlabel('e0')
                    title('Eccentricity - e0')
                    legend('Posterior', 'True Value')

                    subplot(2,3,6)
                    hist(samples(:,6), 30, 'FaceColor', 'b'); hold on
                    line([phi0 phi0], [0 max(histcounts(samples(:,6), 30))./2], 'Color', 'r'); hold on
                    ylabel('P(\phi0|M_0,m_0,z,f_0,e_0)')
                    xlabel('\phi0')
                    title('Phase - \phi0')
                    legend('Posterior', 'True Value')
            end
            function plotOrbit(obj,t,forb,e,logplot)
                %PLOTORBIT Plot orbital evolution parameters
                %
                % Inputs:
                %   t      : Time vector
                %   forb   : Orbital frequency
                %   e      : Eccentricity
                %   logplot: Logarithmic scale flag
            figure
            subplot(1,3,1)
            if logplot==1
                loglog(t,forb,'k');hold on
                set(gca,'fontsize',12);
                set(gca,'LineWidth',1); 
                set(get(gca,'xlabel'),'Fontsize',14);
                set(get(gca,'ylabel'),'Fontsize',14);
                xlabel('Time(s)','FontSize',14,'Fontname','Time new roman'); 
                ylabel('Frequency','FontSize',14,'Fontname','Time new roman'); 
            else
                plot(t,forb,'k');hold on
                set(gca,'fontsize',12);
                set(gca,'LineWidth',1); 
                set(get(gca,'xlabel'),'Fontsize',14);
                set(get(gca,'ylabel'),'Fontsize',14);
                xlabel('Time(s)','FontSize',14,'Fontname','Time'); 
                ylabel('Frequency','FontSize',14,'Fontname','Time new roman'); 
            end

            subplot(1,3,2)
            if logplot==1
                loglog(t,e,'k');hold on
                set(gca,'fontsize',12);
                set(gca,'LineWidth',1); 
                set(get(gca,'xlabel'),'Fontsize',14);
                set(get(gca,'ylabel'),'Fontsize',14);
                xlabel('Time(s)','FontSize',14,'Fontname','Time New roman'); 
                ylabel('Eccentricity','FontSize',14,'Fontname','Time new roman'); 
            else
                plot(t,e,'k');hold on
                set(gca,'fontsize',12);
                set(gca,'LineWidth',1); 
                set(get(gca,'xlabel'),'Fontsize',14);
                set(get(gca,'ylabel'),'Fontsize',14);
                xlabel('Time(s)','FontSize',14,'Fontname','Time'); 
                ylabel('Eccentricity','FontSize',14,'Fontname','Time new roman'); 
            end
             subplot(1,3,3)
            if logplot==1
                loglog(forb,e,'k');hold on
                set(gca,'fontsize',12);
                set(gca,'LineWidth',1); 
                set(get(gca,'xlabel'),'Fontsize',14);
                set(get(gca,'ylabel'),'Fontsize',14);
                xlabel('Frequency','FontSize',14,'Fontname','Time New roman'); 
                ylabel('Eccentricity','FontSize',14,'Fontname','Time new roman'); 
            else
                plot(forb,e,'k');hold on
                set(gca,'fontsize',12);
                set(gca,'LineWidth',1); 
                set(get(gca,'xlabel'),'Fontsize',14);
                set(get(gca,'ylabel'),'Fontsize',14);
                xlabel('Frequency','FontSize',14,'Fontname','Time'); 
                ylabel('Eccentricity','FontSize',14,'Fontname','Time new roman'); 
            end
            end

            function h=tight_subplot(obj,m, n, row, col, gutter, margins, varargin)
                %TIGHT_SUBPLOT Create subplots with minimal whitespace
                %
                % Custom alternative to MATLAB's subplot with better spacing control
                    if nargin<5 || isempty(gutter)
                        gutter = [.002, .002];
                    end
                    if length(gutter)==1
                        gutter(2)=gutter;
                    elseif length(gutter) > 2
                        error('GUTTER must be of length 1 or 2')
                    end          
                    if nargin<6 || isempty(margins)
                        margins = [.06 .01 .04 .04];
                    end         
                    Lmargin = margins(1);
                    Rmargin = margins(2);
                    Bmargin = margins(3);
                    Tmargin = margins(4);          
                    unit_height = (1-Bmargin-Tmargin-(m-1)*gutter(2))/m;
                    height = length(row)*unit_height + (length(row)-1)*gutter(2);            
                    unit_width = (1-Lmargin-Rmargin-(n-1)*gutter(1))/n;
                    width = length(col)*unit_width + (length(col)-1)*gutter(1);         
                    bottom = (m-max(row))*(unit_height+gutter(2))+Bmargin;
                    left   = (min(col)-1)*(unit_width +gutter(1))+Lmargin;              
                    pos_vec= [left bottom width height];        
                    h=subplot('Position', pos_vec, varargin{:});
            end

            function formatted_title = format_title(obj,parameter, median, lower, upper)
                %FORMAT_TITLE Create formatted title with error bars
                %
                % Example output: "Parameter=1.23^{+0.45}_{-0.67}"
                        if isempty(median)
                            median_str = 'N/A';
                        else
                            try
                                median_str = scientific_format2(median)
                            catch
                                median_str = 'N/A';
                            end
                        end
                        if isempty(lower)
                            lower_str = 'N/A';
                        else
                            lower_str = scientific_format(lower);
                        end
                        if isempty(upper)
                            upper_str = 'N/A';
                        else
                            upper_str = scientific_format(upper);
                        end
                        formatted_title = sprintf('%s=%s^{+%s}_{-%s}', parameter, median_str, upper_str, lower_str);
            end
    end
end
