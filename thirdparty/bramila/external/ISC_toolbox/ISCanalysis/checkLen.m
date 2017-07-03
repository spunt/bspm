function checkLen(val)


if val > 1

    n(1)=12;
    for k = 2:15
        n(k) = n(k-1)*2;
    end
    disp(' ')
    disp(['Time-series of length ' num2str(n(val)) ' is required ' ...
        'when the number of frequency subbands is ' num2str(val) '.'])

end