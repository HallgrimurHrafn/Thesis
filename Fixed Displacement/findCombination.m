function valueList=findCombination(target,list,freeOptimize)
% initialize
    valueList=[inf,inf,inf];
%     min 0 spacers used, max 5.
    for n=0:5
%         all possible combinations from list given number of items used.
        combos=combnk(list,n);
        if freeOptimize
            [~,bestInd]=min(abs(sum(combos,2)-target));
            bestCurrentCombo=sum(combos(bestInd,:),2);
        else
    %         targets is upper limit. Thus max(combo*(combo-target<=0)), finds
    %         the closest value to target, excluding larger values.
            [bestCurrentCombo,bestInd]=max(sum(combos,2).*(sum(combos,2)-target<=0));
        end
%         is the newest value better
        if abs(target-bestCurrentCombo)<abs(target-sum(valueList,2))
            valueList=combos(bestInd,:);
        end
%         is the current target good enough? fewer items is better.
        if abs(target-bestCurrentCombo)<=0.015
            break
        end
    end
end