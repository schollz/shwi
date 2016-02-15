function [L0] = getLc2(P,x,F)
    kT=4.1;
    myfun = @(L) F*P/kT-1/(4*(1-x/L).^2)+1/4-x/L;
    L0=fzero(myfun,round(x)*1.5);
end