function conf = plotConfusionMatrix(target, out)
%     conf_mat = confusionmat(target, out);
    t = zeros(length(target),1);
    o = zeros(length(t),1);
    for i = 1:length(t)
        if strcmp(target{i}, 'paper')==1
            t(i) = 0;
        elseif strcmp(target{i}, 'scissor')==1
            t(i) = 1;
        elseif strcmp(target{i}, 'rock')==1
            t(i) = 2;
        end
    end
    
    for i = 1:length(o)
        if strcmp(out{i}, 'paper')==1
            o(i) = 0;
        elseif strcmp(out{i}, 'scissor')==1
            o(i) = 1;
        elseif strcmp(out{i}, 'rock')==1
            o(i) = 2;
        end
    end
    
    conf = confusionmat(t, o);
end