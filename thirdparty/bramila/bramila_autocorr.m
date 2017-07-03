function df=bramila_autocorr(x,y)
    % See appendix B in http://www.sciencedirect.com/science/article/pii/S1053811911013000
    % Usage df=bramila_autocorr(x,y)
    %   x = time series (column vector)
    %   y = second time series (column vector)
    %   df = estimated degrees of freedom
    % For autocorrelation just run df=bramila_autocorr(x,x)

    N=length(x);
    j=round(N/5);
    a=autocorr(x,j);
    b=autocorr(y,j);
    temp=0;
    for lag=1:j
        temp=temp+(2/N)*((N-lag)/N)*a(lag+1)*b(lag+1);
    end
    temp=1/N+temp;
    df=1/temp;
    df=df;
    
