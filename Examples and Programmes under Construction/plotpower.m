function [pwWFSim,timeWFSim,pwLES,timeLES] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, plotLES );

sol          = sol_array(end);
measuredData = sol.measuredData;

solArrayTurb = [sol_array.turbine];
pwWFSim      = [solArrayTurb.power];
timeWFSim    = Wp.sim.time(2:sol.k+1);
timeLES      = Wp.sim.time(1:Wp.sim.NN);

pwLES = [];
if plotLES == 1
    for i = 1:Wp.sim.NN
        pwLES(:,i) = sol_array(i).measuredData.power;
    end
end

if scriptOptions.powerForecast > 0
    k_end = min(Wp.sim.NN,scriptOptions.powerForecast+sol.k);
    disp(['Forecasting power for k = ' num2str(sol.k+1) ' until k = ' num2str(k_end) '.']);
    sol_tmp    = sol;  % Copy current sol
    sol_tmp.uu = sol.u;
    sol_tmp.vv = sol.v;
    sol_tmp.pp = sol.p;
    timeFC = []; pwFC = [];
    while sol_tmp.k < k_end
        [ sol_tmp,~ ] = WFSim_timestepping( sol_tmp, sys, Wp, scriptOptions );
        pwFC(:,sol_tmp.k-sol.k) = sol_tmp.turbine.power; % Extract power
        timeFC(sol_tmp.k-sol.k) = sol_tmp.time;
    end
    pwFC = movmean(pwFC',[10 0])'; % Low-pass filter power forecast
    t_st = min([scriptOptions.powerForecast, 60]); % t-shortterm
    try
        RMSE_FC_shortterm = sqrt(mean((pwLES(:,timeFC(1:t_st))-pwFC(:,1:t_st)).^2,2));
        RMSE_FC_longterm  = sqrt(mean((pwLES(:,timeFC)-pwFC).^2,2));
        %             RMSE_FC_shortterm = mean(abs(pwLES(:,timeFC(1:t_st))-pwFC(:,1:t_st)),2);
        %             RMSE_FC_longterm  = mean(abs(pwLES(:,timeFC)-pwFC),2);
        disp(['Short-term Mean RMSE forecasted vs. true power for k = ' num2str(sol.k+1) ':' num2str(sol.k+t_st) ' is ' num2str(mean(RMSE_FC_shortterm),'%10.2e\n') '.']);
        disp(['Long-term  Mean RMSE forecasted vs. true power for k = ' num2str(sol.k+1) ':' num2str(k_end) ' is ' num2str(mean(RMSE_FC_longterm),'%10.2e\n') '.']);
    end
end

subplotDim = numSubplots(Wp.turbine.N); % Determine optimal layout for subplots
for j = 1:Wp.turbine.N
    subplot(subplotDim(1),subplotDim(2),j);
    if plotLES == 1
        plot(timeLES,pwLES(j,:),'k-','lineWidth',0.25,'DisplayName', ['LES']); hold on
    end
    plot(timeWFSim,pwWFSim(j,:),'-','lineWidth',1.0,'DisplayName', ['WFObs: ' strucObs.filtertype]); hold on;
    if scriptOptions.powerForecast > 0 && length(timeFC) > 0
    plot(timeFC,pwFC(j,:),'-','lineWidth',0.75,'DisplayName',['WFObs: ' strucObs.filtertype ' (FC)']); hold on;
    end
    legend('-DynamicLegend');
    xlabel('Time (s)');
    xlim([0 Wp.sim.NN]);
    if plotLES == 1
        ylim([0 1e6*ceil(max([pwWFSim(:); pwLES(:)])/1e6)]);
    else
        ylim([0 1e6*ceil(max([pwWFSim(:)])/1e6)]);
    end
    grid on;
    ylabel('Power (W)');
    title(['Turbine ' num2str(j) '']);
end;