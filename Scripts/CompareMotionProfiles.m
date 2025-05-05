function CompareMotionProfiles(profiles, labels)
% COMPAREMOTIONPROFILES – confronta più motion profiles
%
%   compareMotionProfiles(profiles, labels)
%
% Input:
%   profiles: cell array di struct con campi time, position, velocity, acceleration
%   labels : cell array di stringhe per la legenda
%
% Esempio di utilizzo: 

% costruisco i profili q2 e q3
% [t2, s2, v2, a2] = restToRestMotion(L2, A2, V2, 1);
% [t3, s3, v3, a3] = restToRestMotion(L3, A3, V3, 1);
% 
% % preparo il cell array
% profiles = {
%   struct('time',t2,'position',s2,'velocity',v2,'acceleration',a2), ...
%   struct('time',t3,'position',s3,'velocity',v3,'acceleration',a3)
% };
% 
% labels = {'q2','q3'};
% 
% % comparo
% compareMotionProfiles(profiles, labels);


  n = numel(profiles);
  colors = lines(n);

  figure('Name','Confronto Motion Profiles','Units','normalized','Position',[.2 .3 .5 .5]);

  % Plot s(t)
  subplot(3,1,1); hold on; grid on;
  for i=1:n
    plot(profiles{i}.time, profiles{i}.position, 'Color', colors(i,:), 'LineWidth',1.5);
  end
  title('s(t)'); xlabel('t [s]'); ylabel('s(t)');
  legend(labels,'Location','best');

  % Plot v(t)
  subplot(3,1,2); hold on; grid on;
  for i=1:n
    plot(profiles{i}.time, profiles{i}.velocity, 'Color', colors(i,:), 'LineWidth',1.5);
  end
  title('v(t)'); xlabel('t [s]'); ylabel('v(t)');
  legend(labels,'Location','best');

  % Plot a(t)
  subplot(3,1,3); hold on; grid on;
  for i=1:n
    plot(profiles{i}.time, profiles{i}.acceleration, 'Color', colors(i,:), 'LineWidth',1.5);
  end
  title('a(t)'); xlabel('t [s]'); ylabel('a(t)');
  legend(labels,'Location','best');
end