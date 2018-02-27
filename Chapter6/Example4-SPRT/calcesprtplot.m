function calcesprtplot(testdata, alarm, islabel)

if strcmp(islabel, 'true')
    label = testdata(:,end);
    falsealarmposcandidates = find(alarm==1);
    if isempty(falsealarmposcandidates)
        falsealarmpos = [];
    end 
    missedalarmposcandidates = find(alarm==0);
    if isempty(missedalarmposcandidates)
        missedalarmpos = [];
    end
    falsealarmpospt = 1;
    missedalarmpospt = 1;
    for m=1 : length(falsealarmposcandidates)
        pos = falsealarmposcandidates(m);
        if label(pos)==0
            falsealarmpos(falsealarmpospt) = pos;
            falsealarmpospt = falsealarmpospt+1;
        end
        if ~exist('falsealarmpos')
            falsealarmpos = [];
        end
    end
    for m=1 : length(missedalarmposcandidates)
        pos = missedalarmposcandidates(m);
        if label(pos)==1
            missedalarmpos(missedalarmpospt) = pos;
            missedalarmpospt = missedalarmpospt+1;
        end
        if ~exist('missedalarmpos')
            missedalarmpos = [];
        end
    end 

    figure,
    plot(find(label==0), testdata(find(label==0)), 'bo');
    hold on, plot(find(label==1), testdata(find(label==1)), 'ro');
    hold on, plot(falsealarmpos, testdata(falsealarmpos), 'go');
    hold on, plot(missedalarmpos, testdata(missedalarmpos), 'ko');
    xlabel('# of data points'), ylabel('Magnitude');
    if isempty(falsealarmpos) & isempty(missedalarmpos)
        legend('Healthy data points', 'Anomalous data points');
    elseif isempty(falsealarmpos) & ~isempty(missedalarmpos)
        missedalarmrate = (length(missedalarmpos)/length(find(label==1)))*100;
        str = {'Missed alarm rate = ', [num2str(missedalarmrate) ' %']};
        mxpos = ceil(median(find(label==1)));
        mypos = ceil(median(testdata(find(label==1))));
        text(mxpos, mypos, str, 'Color', 'red', 'FontSize', 14);
        legend('Healthy data points', 'Anomalous data points',...
            'Missed alarmed data points');
    elseif ~isempty(falsealarmpos) & isempty(missedalarmpos)
        falsealarmrate = (length(falsealarmpos)/length(find(label==0)))*100;
        str = {'False alarm rate = ', [num2str(falsealarmrate) ' %']};
        fxpos = ceil(median(find(label==0)));
        fypos = ceil(median(testdata(find(label==0))));
        text(fxpos, fypos, str, 'Color', 'red', 'FontSize', 14);
        legend('Healthy data points', 'Anomalous data points',...
            'False alarmed data points');
    elseif ~isempty(falsealarmpos) & ~isempty(missedalarmpos)
        falsealarmrate = (length(falsealarmpos)/length(find(label==0)))*100;
        missedalarmrate = (length(missedalarmpos)/length(find(label==1)))*100;
        fstr = {'False alarm rate = ', [num2str(falsealarmrate) ' %']};
        mstr = {'Missed alarm rate = ', [num2str(missedalarmrate) ' %']};
        mxpos = ceil(median(find(label==1)));
        mypos = ceil(median(testdata(find(label==1))));
        fxpos = ceil(median(find(label==0)));
        fypos = ceil(median(testdata(find(label==0))));
        text(mxpos, mypos, mstr, 'Color', 'red', 'FontSize', 14);
        text(fxpos, fypos, fstr, 'Color', 'red', 'FontSize', 14);
        legend('Healthy data points', 'Anomalous data points',...
            'False alarmed data points', 'Missed alarmed data points');
    end
    
elseif ~strcmp(islabel, 'true');
    alarmpos = find(alarm==1);
    figure,
    plot(testdata, 'bo');
    hold on, plot(alarmpos, testdata(alarmpos), 'ro');
    xlabel('# of data points'), ylabel('Magnitude');
    legend('Test data points', 'Anomalous data points');
end

end

