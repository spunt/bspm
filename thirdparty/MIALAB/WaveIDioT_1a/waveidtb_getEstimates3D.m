function [y_est] = waveidtb_getEstimates3D(D , Dprev , sig , K)

% Initialize the classification
% x = zeros( size(D) );
y = Dprev;

% Calculate the mask
x = (abs(D).*abs(y) ) >= ( (K^2)*sig);

% Estimate the wavelet coefficients after the mask has been obtained
% Get Magnitudes
m = abs(D);

% Compute the locally averaged coefficient magnitude
[M,N,O] = size(D);

% Find the ratio of 1's to 0's in x
if ( sum(x(:)) < M*N*O )
    r = sum(x(:))/(size(x,1) * size(x,2) * size(x,3) - sum(x(:)));
else
    r = 100;
end;

if (r < 1e-4)
    
    y_est = zeros(M,N);
    
elseif ( max( m(x==0) ) < 1 )
    
    y_est = D;
    
else
    
    % Calcualte the Mximum Likelihood histograms
    % Smart way to compute 3-D indices of locations of 1s and 0s without
    % ind2sub
    %     [xS0,yS0] = find(reshape(x , M,N*O)==0);
    %     zS0 = ceil(yS0/N);
    %     [xS1,yS1] = find(reshape(x , M,N*O)==1);
    %     zS1 = ceil(yS1/N);
    
    ind0 = find(x==0);
    [xS0,yS0,zS0] = ind2sub(size(x),ind0);
    ind1 = find(x==1);
    [xS1,yS1,zS1] = ind2sub(size(x),ind1);
    
    % Find the probability maps and Rearrange probabilities to maps;
    
    % Compute p( x | 1 )
    %     bins = ceil(max(abs(D(:))));
        P_m_x_1 =  hist( m(x==1) , numel( m( x==1 ) ) ) ./numel(m);
%     m_x_1 = m(:).*(x(:)>0);
%     P_m_x_1 = ( hist( nonzeros(m_x_1) , numel( nonzeros(m_x_1) ) ) )./numel(m);
    P_m_x_1Map = zeros( size(m) );
    P_m_x_1Map(ind1) = P_m_x_1;
    
    
    % Compute p( x | 0 )
        P_m_x_0 = ( hist( m( x==0 ) , numel( m( x == 0 ) ) ) )./ numel(m);
%     P_m_x_0 = ( hist( nonzeros(m(:).*(x(:)==0)) , numel( nonzeros(m.*(x==0)) ) ) )./numel(m);
    P_m_x_0Map = zeros( size(m) );
    P_m_x_0Map(ind0) = P_m_x_0;
    
    
    % Check the number of elements that are classified as noisy
    if ( numel(nonzeros(P_m_x_0Map)) < 5 )
        y_est = D;
    else
        
        % Calculate the apriori maps using Coarser details
        % Trick used by replacing the D1 with D2 where M1 = 1 and scale
        
        % Find the apriori probabilities
        
        % Signal
        
        D1 = waveidtb_neighborhood3d( padarray( D , [1 1 1] , 'both' ) , 27);
%         D1 = waveidtb_neighborhood3d_mex( padarray( D , [1 1 1] , 'both' ));
        
        D2 = waveidtb_neighborhood3d( padarray( Dprev , [1 1 1] , 'both' ) , 27);
%         D2 = waveidtb_neighborhood3d_mex( padarray( Dprev , [1 1 1] , 'both' ));
        
        % Note the difference here
        temp1 = waveidtb_neighborhood3d( padarray( x , [1 1 1] , 'both' ) , 27);
%         temp1 = waveidtb_neighborhood3d_mex( padarray( double(x) , [1 1 1] , 'both' ));
        
        temp2 = repmat( temp1(14,:) , [27 1] );  % P(Xk==1)
        clear temp1;
        D1(14,:) = D2(14,:);      % Replace the D with Dprev element
        e = reshape ( sqrt( sum( (temp2.*D1 ).^2 , 1) ./ 27 ) , size(D)+2 );
        e = e(2:end-1,2:end-1,2:end-1);
        
        P_e_x_1 = ( hist( e( x==1 ) , numel( e(x==1) ) ) )./numel(e);
%         P_e_x_1 = ( hist( nonzeros(e(:).*(x(:)==1)) , numel( nonzeros(e.*(x==1)) ) ) )./numel(e);
        P_e_x_1Map = zeros( size( e ) );
        P_e_x_1Map(ind1) = P_e_x_1;
        clear temp2;
        
        %Noise
        
        temp1 = waveidtb_neighborhood3d( padarray( 1 - x , [1 1 1] , 'both' ) , 27);
%         temp1 = waveidtb_neighborhood3d_mex( padarray( double(1 - x) , [1 1 1] , 'both' ));
        
        temp2 = repmat( temp1(14,:) , [27 1] );  % P(Xk==0)
        clear temp1;
        D1(14,:) = D2(14,:);      % Replace the D with Dprev element
        e = reshape ( sqrt( sum( (temp2.*D1 ).^2 , 1) ./ 27 ) , size(D)+2 );
        e = e(2:end-1,2:end-1,2:end-1);
        
        P_e_x_0 = ( hist( e(x==0 ) , numel( e(x==0) ) ) )./numel(e);
%         P_e_x_0 = ( hist( nonzeros(e(:).*(x(:)==0)) , numel( nonzeros(e.*(x==0)) ) ) )./numel(e);
        P_e_x_0Map = zeros( size( e ) );
        P_e_x_0Map(ind0) = P_e_x_0;
        
        % Perform the Shrinkage Procedure
        
        % Poisson Distribution ratio maps (Seen from the distributions of
        % signal and noise)
        P_m_x_0Map(P_m_x_0Map==0) = min(nonzeros(P_m_x_0Map(:)));
        P_e_x_0Map(P_e_x_0Map==0) = min(nonzeros(P_m_x_1Map(:)));
        xeta = (P_m_x_1Map)./(P_m_x_0Map);      % Likelihood ratio
        neta = (P_e_x_1Map)./(P_e_x_0Map);      % apriori ratio
        factor = (r.*xeta.*neta);
        %         *10^round(abs(log10(max(xeta(:).*neta(:)))));
        y_est = D.*(factor./(1+factor));
        y_est(isnan(y_est)) = 0;
    end;
    
end;
