function seq = introduce_injection_delay(seq, aif_mask, injection_delay)
% Assuming the input DSA series is 3D, where the first dimension is time.
% The value range is [0, 255], and pixel values are higher with more
% contrast.

    aif=mean(seq(:,aif_mask==1),2);
    n=numel(aif);

    % --- gaussian smoothing ---
    w = gausswin(15);
    w = w/sum(w);
    aif = filter(w, 1, aif);
    % ----------------------------
    
    [~, uptake_idx] = max(diff(aif) + circshift(diff(aif), -1) + circshift(diff(aif), -2));
    if injection_delay < uptake_idx
        seq = seq((uptake_idx-injection_delay):end, :, :);
    elseif injection_delay > uptake_idx
        seq = padarray(seq, [injection_delay-uptake_idx, 0, 0], 'replicate', 'pre');
    end
end