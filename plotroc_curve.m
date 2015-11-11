function [fpr1 tpr1 ] = plotroc_curve( input_args )
%PLOTTROC_CURVE Summary of this function goes here
%   Detailed explanation goes here
score = input_args(:,1);
label = input_args(:,2);

maxscore = max(score);
minscore = min(score);

if maxscore > minscore
    score = (score - minscore)/maxscore;
else
    score = 0*score;
    maxscore = 0.;
    minscore = score(1);
end

% Expand score and label
targets = [label'; 1-label'];
outputs = [1-score'; score'];
[tpr,fpr,thresholds] = roc(targets,outputs);

fig = figure(1);
plot(fpr{1}, tpr{1}, 'g-');
fpr1 = fpr{1};
tpr1 = tpr{1};
hold on;
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC');
grid on;
hold off;
dcm_obj = datacursormode(fig);
myupdatefcn1 = @(obj,event_obj)myupdatefcn(obj,event_obj,(1-thresholds{1})*maxscore+minscore, fpr{1}); 
set(dcm_obj,'UpdateFcn',myupdatefcn1);
hold on
plot([0 1],[1 0], 'b-.');
hold off
end



function output_txt = myupdatefcn(obj,event_obj,thresholds, fpr)
% obj              Currently not used (empty)
% event_obj        Handle to event object
% output_txt       Data cursor text string (string or cell array of
%                  strings)
pos = get(event_obj,'Position');
t = 0;
for i=1:length(fpr)
    if abs(fpr(i) - pos(1)) <1.e-6
        t = thresholds(i);
        break;
    end
end
output_txt =  {['fpr: ',num2str(pos(1))],...
       ['tpr: ',num2str(pos(2))],...
       ['cut: ',num2str(t) ]};
end