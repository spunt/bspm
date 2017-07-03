function data = bramila_derivative(data,order)
    % INPUT:
    %       data = separate time-series in columns (important!)
    %       order = derivative order
    % OUTPUT:
    %       res = derivative of proper order
    order = round(order);    

%     data = data'; % gradient operates on rows
%     for i=1:order % repeat until wanted order
%         data = gradient(data);
%     end
%     data = data';

    
    dts=diff(data,order);
    data=[zeros(order,size(dts,2)); dts];  % first element is a zero, as per Power et al 2014
    
    
    % edges (first and last element) are approximated by forward/backward differences
    % all other elements by centered differences

end