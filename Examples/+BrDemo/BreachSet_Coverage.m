% Demonstrates basic operations relative to coverage estimates

%%
% Creates some set with some 9 arbitrary parameters

params = {'p1', 'p2', 'p3', ...
          'q1', 'q2', 'q3', ...
          'r1', 'r2', 'r3', ...
         };
B0 = BreachSet(params);

%% 
% Set unit range for all parameters
B3 = B0.copy();
B3.SetParamRanges(params(1:3), [0 1]);
B6 = B0.copy();
B6.SetParamRanges(params(1:6), [0 1]);
B9 = B0.copy();
B9.SetParamRanges(params, [0 1]);

il = 100:100:1000;
epsi = 0.1;
delta = 0.2; 

%% 
% 
profile on
[covocc3, covlog3, covent3 ] = compute_batch_coverage(B3, il, epsi, delta); 
profile report
%%
%
[covocc6, covlog6, covent6 ] = compute_batch_coverage(B6, il, epsi, delta); 

%%
%
[covocc9, covlog9, covent9 ] = compute_batch_coverage(B9, il, epsi, delta); 

%%
%
figure;
grid on;
hold on;
% plot(il, covocc3, 'b', il, covlog3,'+-b',il, covent3,'*-b',...
%      il, covocc6,'g' ,il, covlog6,'+-g',il, covent6,'*-g',...
%      il, covocc9,'r', il, covlog9,'+-r',il, covent9,'*-r');
plot(il, covocc3, 'b', il, covlog3,'+-b', ...
     il, covocc6,'g' ,il, covlog6,'+-g', ...
     il, covocc9,'r', il, covlog9,'+-r');


% legend('Occ(3d)','LogOcc(3d)','Ent(3d)',...
%     'Occ(6d)','LogOcc(6d)','Ent(6d)',...
%     'Occ(9d)','LogOcc(9d)','Ent(6d)');
legend('Occ(3d)','LogOcc(3d)',...
    'Occ(6d)','LogOcc(6d)',...
    'Occ(9d)','LogOcc(9d)');

%%
% Quasi-random vs Uniform-Random

%% 
% 

[covocc, covlog] = compute_batch_coverage(B9, il, epsi, delta); 






