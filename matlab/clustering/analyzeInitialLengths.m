function analyzeInitialLengths(r)
    %% Analyze initial length2
    l1s = [];
    peaks = [];
    allInitialLengths = [];
    allInitialLengthsError = [];
    allNumPeaks =[];
    for i=1:length(r)
        peaknum =  length(r{i}.F)-1;
        l1 = r{i}.L(2)-r{i}.L(1);
        if (l1>0 && l1<150)   
            peaks = [peaks;  peaknum];
            l1s = [l1s; l1];
        end
    end
    plot(peaks,l1s,'o')
    for peaknum=3:8
        for foocount = 1:100
            [mu_e,sigma_e,p_e]=gaussian_mixture_model(l1s(peaks==peaknum)',3);
            [m,mind] = max(p_e);
            mu_best = mu_e(mind);
            sigma_best =  sigma_e(mind);
            if isnan(mu_best) && isnan(sigma_best)
                [mu_e,sigma_e,p_e]=gaussian_mixture_model(l1s(peaks==peaknum)',2);
                [m,mind] = max(p_e);
                mu_best = mu_e(mind);
                sigma_best =  sigma_e(mind);
            end
            if ~isnan(mu_best) && ~isnan(sigma_best)
                goodind = find(p_e>0.1);
                allInitialLengths = [allInitialLengths; mu_e(goodind);];
                allInitialLengthsError = [allInitialLengthsError; sigma_e(goodind);];
                allNumPeaks = [allNumPeaks; ones(length(goodind),1)*peaknum];
                break
            end
        end
    end
    figure(13)
    hold on;
    plot(allNumPeaks,allInitialLengths,'or')
end
