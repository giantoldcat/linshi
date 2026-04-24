% robot3Trajectory_simulink.m
%% MATLAB 高保真数字孪生仿真：终极完整版（去掉直线IK，保留工作区导出）
clear; clc; close all;

%% ===================== 可复现设置 =====================
rng(20260408,'twister');

%% ===================== 日志初始化 =====================
BO_LOG = struct('iter',[],'f',[],'L',[],'Lref',[],'w',[],'round',[]);
setappdata(0,'BO_LOG',BO_LOG);
setappdata(0,'BO_ROUND_ID',0);
setappdata(0,'BO_TMP',struct('L',[],'Lref',[]));

%% ===================== 0. 数字孪生参数 =====================
twin.online_mode = 'safety_balanced';
twin.safety_margin = 0.04;
twin.friction_scale_j = ones(1,6);
twin.friction_min = 0.85;
twin.friction_max = 1.20;

twin.ground_clearance = 0.05;
twin.ground_penalty_k = 2e4;

twin.joint_limit = pi;
twin.joint_limit_penalty_k = 1e4;

plant.bias = [0.010, -0.008, 0.012, 0.005, -0.006, 0.004];
plant.drift_rate = [0.002, -0.001, 0.0015, 0.001, -0.0012, 0.0008];
plant.noise_std = 0.004;
plant.delay_steps = 2;

lambda_fx = 0.60;

%% ===================== 1. 机械臂模型 =====================
p = [
    -0.019225  -0.000523   0.100684;
    -0.0155     0.035      0.0285;
     0.03665    0          0.146;
    -0.01855   -0.01       0.052;
     0.0168     0.127      0;
    -0.018506   0.077     -0.000106
];
ax = [
     0  0 -1;
    -1  0  0;
    -1  0  0;
     0 -1  0;
    -1  0  0;
     0 -1  0
];

robot = rigidBodyTree('DataFormat','row','MaxNumBodies',6);
parentName = robot.BaseName;
for i = 1:6
    body = rigidBody(sprintf('link%d', i));
    jnt  = rigidBodyJoint(sprintf('Joint%d', i), 'revolute');
    jnt.JointAxis = ax(i,:);
    setFixedTransform(jnt, trvec2tform(p(i,:)));
    body.Joint = jnt;
    addBody(robot, body, parentName);
    parentName = body.Name;
end

%% ===================== 2. 起点 + 笛卡尔目标点IK =====================
q_start = [0,0,0,0,0,0];
p_target = [0.10, 0.22, 0.28];

T0 = getTransform(robot, q_start, 'link6');
R_target = T0(1:3,1:3);

T_target = eye(4);
T_target(1:3,1:3) = R_target;
T_target(1:3,4) = p_target(:);

[q_end, okIK, infoIK] = ik_solve_with_rbt(robot, q_start, T_target);

if ~okIK
    searchOpt.range = 0.03;
    searchOpt.step = 0.01;
    searchOpt.keepSameZFirst = true;

    [q_end, p_target_used, okSearch, infoSearch] = ...
        ik_fallback_nearby_search(robot, q_start, p_target, R_target, searchOpt);

    if ~okSearch
        error('IK失败且邻域±3cm未找到可达点。请调整目标点。');
    else
        disp('原目标IK失败，已切换到邻域最近可达点。');
        disp(['原目标 = [', num2str(p_target,'%.4f '), ']']);
        disp(['替代目标 = [', num2str(p_target_used,'%.4f '), ']']);
        disp(['|dp| = ', num2str(norm(p_target_used-p_target),'%.4f'), ' m']);
        disp(['IK误差 = ', num2str(infoSearch.bestPosErr,'%.6f'), ' m']);
    end
else
    p_target_used = p_target;
    disp(['IK成功，目标点 = [', num2str(p_target_used,'%.4f '), ']']);
    disp(['IK误差 = ', num2str(infoIK.posErr,'%.6f'), ' m']);
end

q_end = max(-pi, min(pi, q_end));
disp(['q_end(deg) = [', num2str(rad2deg(q_end),'%.2f '), ']']);

%% ===================== 3. 障碍物 =====================
q_mid_linear = (q_start + q_end)/2;
T_mid = kinematics_URDF_model(q_mid_linear, p, ax);
base_obs_center = T_mid(1:3,4)';

obs.centers = [
    base_obs_center;
    0.06,  -0.01, 0.20;
    0.003,  0.25, 0.20
];
obs.radii = [0.035; 0.035; 0.045];
obs.num = size(obs.centers,1);

%% ===================== 4. 优化设置 =====================
t_total = 5;
q_via_guess = q_mid_linear;

vars = [
    optimizableVariable('x1',[-5,-1],'Type','real')
    optimizableVariable('x2',[-1,1.4],'Type','real')
    optimizableVariable('x3',[-4,0],'Type','real')
];

inner_opt = optimset('Display','off','MaxIter',140,'MaxFunEvals',520,'TolX',8e-3,'TolFun',5e-4);

%% ===================== 5. 第1轮 =====================
disp('================================================================');
disp('第1轮：BayesOpt + 内层优化');
setappdata(0,'BO_ROUND_ID',1);
setappdata(0,'BO_TMP',struct('L',[],'Lref',[]));

objFun1 = @(X) bo_outer_objective_twin(X, q_start, q_end, t_total, p, ax, obs, q_via_guess, twin);
try
    bo_results1 = bayesopt(objFun1, vars, ...
        'MaxObjectiveEvaluations', 15, ...
        'IsObjectiveDeterministic', true, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'Verbose', 0, 'PlotFcn', {}, 'OutputFcn', @bo_output_logger);

    xb1 = bo_results1.XAtMinObjective;
    w_best1 = [10^(xb1.x1), 10^(xb1.x2), 10^(xb1.x3)];
catch ME
    warning('BO:Round1Failed','%s',ME.message);
    w_best1 = [5e-3, 2.0, 2e-1];
end

[q_via_opt1, fval1, exitflag1, out1] = fminsearch(@(qv) cost_function_weighted_twin(...
    qv, q_start, q_end, t_total, p, ax, obs, w_best1, twin), q_via_guess, inner_opt);

if exitflag1 <= 0
    warning('OPT:Round1InnerNoConverge','第1轮内层未严格收敛，f=%.6f',fval1);
end

disp('---------------- 第1轮内层优化信息 ----------------');
disp(['w_best1 = [', num2str(w_best1,'%.3e '), ']']);
disp(['exitflag = ', num2str(exitflag1), ', iterations = ', num2str(out1.iterations), ...
      ', funcCount = ', num2str(out1.funcCount), ', fval = ', num2str(fval1,'%.6f')]);
disp('---------------------------------------------------');

[q_opt_t1, dq_opt_t1, ~] = trajectory_SDPOA(q_start, q_via_opt1, q_end, t_total, 0:1e-4:t_total);

%% ===================== 6. 同步校准 =====================
disp('开始实体数据同步与参数校准...');
[q_meas_t, dq_meas_t] = simulate_real_plant(q_opt_t1, dq_opt_t1, plant, 1e-4);
err = sync_real_data(q_opt_t1, dq_opt_t1, q_meas_t, dq_meas_t);
twin = update_twin_params(twin, err);

disp(['校准后 fixed safety_margin = ', num2str(twin.safety_margin,'%.5f'), ...
      ', friction_scale_j = [', num2str(twin.friction_scale_j,'%.3f '), ']']);

%% ===================== 7. 第2轮（按日志f-L权衡） =====================
disp('第2轮：BayesOpt + 按日志f-L权衡选解');
setappdata(0,'BO_ROUND_ID',2);
setappdata(0,'BO_TMP',struct('L',[],'Lref',[]));

objFun2 = @(X) bo_outer_objective_twin(X, q_start, q_end, t_total, p, ax, obs, q_via_opt1, twin);
try
    bo_results2 = bayesopt(objFun2, vars, ...
        'MaxObjectiveEvaluations', 15, ...
        'IsObjectiveDeterministic', true, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'Verbose', 0, 'PlotFcn', {}, 'OutputFcn', @bo_output_logger);

    X2 = bo_results2.XTrace;
    BO_LOG = getappdata(0,'BO_LOG');
    idxR2 = find(BO_LOG.round == 2);

    nUse = min(height(X2), numel(idxR2));
    fLog2 = BO_LOG.f(idxR2(1:nUse));
    LLog2 = BO_LOG.L(idxR2(1:nUse));

    idxPool = (1:nUse)';
    fn = normalize01(fLog2(idxPool));
    Ln = normalize01(LLog2(idxPool));
    score = lambda_fx*fn + (1-lambda_fx)*Ln;
    [~,kk] = min(score);
    idxBest2 = idxPool(kk);

    w_best2 = [10^(X2.x1(idxBest2)), 10^(X2.x2(idxBest2)), 10^(X2.x3(idxBest2))];
    q_via_seed2 = q_via_opt1;

    disp('---------------- 第2轮候选权衡结果(与日志一致) ----------------');
    disp(['lambda_fx = ', num2str(lambda_fx,'%.2f')]);
    disp(['选中第2轮 Iter = ', num2str(idxBest2)]);
    disp(['f_log = ', num2str(fLog2(idxBest2),'%.6f'), ', L_log = ', num2str(LLog2(idxBest2),'%.6f')]);
    disp(['w_best2 = [', num2str(w_best2,'%.6e '), ']']);
    disp('----------------------------------------------------------------');
catch ME
    warning('BO:Round2Failed','%s',ME.message);
    w_best2 = w_best1;
    q_via_seed2 = q_via_opt1;
end

[q_via_opt2, fval2, exitflag2, out2] = fminsearch(@(qv) cost_function_weighted_twin(...
    qv, q_start, q_end, t_total, p, ax, obs, w_best2, twin), q_via_seed2, inner_opt);

if exitflag2 <= 0
    warning('OPT:Round2InnerNoConverge','第2轮内层未严格收敛，f=%.6f',fval2);
end

disp('---------------- 第2轮内层优化信息 ----------------');
disp(['w_best2 = [', num2str(w_best2,'%.3e '), ']']);
disp(['exitflag = ', num2str(exitflag2), ', iterations = ', num2str(out2.iterations), ...
      ', funcCount = ', num2str(out2.funcCount), ', fval = ', num2str(fval2,'%.6f')]);
disp('---------------------------------------------------');

%% ===================== 8. 最终轨迹 =====================
q_via_final = q_via_opt2;
t_sample = 1e-4;
t = 0:t_sample:t_total;
N = length(t);

[q_actual_t, dq_actual_t, ddq_actual_t] = trajectory_SDPOA(q_start, q_via_final, q_end, t_total, t);
[q_bad_t, ~, ~] = trajectory_SDPOA(q_start, q_mid_linear, q_end, t_total, t);

%% ===================== 9. 运动学解算 =====================
ds_factor = 50;
ds_idx = 1:ds_factor:N;
ds_N = length(ds_idx);
ds_pos_opt = zeros(ds_N,3);
ds_pos_bad = zeros(ds_N,3);

for i = 1:ds_N
    k = ds_idx(i);
    Ttmp1 = kinematics_URDF_model(q_actual_t(k,:), p, ax);
    ds_pos_opt(i,:) = Ttmp1(1:3,4)';
    Ttmp2 = kinematics_URDF_model(q_bad_t(k,:), p, ax);
    ds_pos_bad(i,:) = Ttmp2(1:3,4)';
end

ds_t = t(ds_idx);
pos_actual_t = zeros(N,3);
pos_actual_t(:,1) = interp1(ds_t, ds_pos_opt(:,1), t, 'spline');
pos_actual_t(:,2) = interp1(ds_t, ds_pos_opt(:,2), t, 'spline');
pos_actual_t(:,3) = interp1(ds_t, ds_pos_opt(:,3), t, 'spline');

pos_bad_t = zeros(N,3);
pos_bad_t(:,1) = interp1(ds_t, ds_pos_bad(:,1), t, 'spline');
pos_bad_t(:,2) = interp1(ds_t, ds_pos_bad(:,2), t, 'spline');
pos_bad_t(:,3) = interp1(ds_t, ds_pos_bad(:,3), t, 'spline');

%% ===================== 9.1 导出6关节优化轨迹到工作区（可直接给Simulink） =====================
traj_t = t(:);
q_opt  = q_actual_t;
dq_opt = dq_actual_t;
ddq_opt= ddq_actual_t;

traj_opt.time = traj_t;
traj_opt.q    = q_opt;
traj_opt.dq   = dq_opt;
traj_opt.ddq  = ddq_opt;
assignin('base','traj_opt',traj_opt);

for j = 1:6
    assignin('base', sprintf('pos_ref%d', j), [traj_t, q_opt(:,j)]);
    assignin('base', sprintf('vel_ref%d', j), [traj_t, dq_opt(:,j)]);
    assignin('base', sprintf('acc_ref%d', j), [traj_t, ddq_opt(:,j)]);

    assignin('base', sprintf('q%d_opt', j),   q_opt(:,j));
    assignin('base', sprintf('dq%d_opt', j),  dq_opt(:,j));
    assignin('base', sprintf('ddq%d_opt', j), ddq_opt(:,j));
end

assignin('base','traj_time',traj_t);
assignin('base','q_opt_all',q_opt);
assignin('base','dq_opt_all',dq_opt);
assignin('base','ddq_opt_all',ddq_opt);

disp('已导出到工作区：pos_ref1~6, vel_ref1~6, acc_ref1~6, traj_opt 等。');

%% ===================== 导出 pos_ref 到 Excel =====================
% 第一列为时间 traj_t，第二到第七列为 1-6 关节的位置数据 q_opt
export_data = [traj_t, q_opt(:,1), q_opt(:,2), q_opt(:,3), q_opt(:,4), q_opt(:,5), q_opt(:,6)];

% 定义 Excel 表头列名
colNames = {'Time', 'pos_ref1', 'pos_ref2', 'pos_ref3', 'pos_ref4', 'pos_ref5', 'pos_ref6'};

% 将数据转换为 Table 格式，以便带有表头写入
export_table = array2table(export_data, 'VariableNames', colNames);

% 导出为 Excel 表格
excel_filename = 'Robot_Trajectory_pos_ref.xlsx';
writetable(export_table, excel_filename);

disp(['已成功将 pos_ref1~6 导出到 Excel 表格：', excel_filename]);

%% ===================== 10. BO日志输出 =====================
BO_LOG = getappdata(0,'BO_LOG');
if ~isempty(BO_LOG.iter)
    disp(' ');
    disp('================== BO 每次更迭日志（两轮） ==================');
    fprintf('Round\tIter\tf(X)\t\t\tL_withObs(m)\tL_ref_noObs(m)\t[w_acc, w_obs, w_len]\n');
    for k = 1:length(BO_LOG.iter)
        fprintf('%d\t%d\t%.6e\t%.6f\t%.6f\t[%.2e, %.2e, %.2e]\n', ...
            BO_LOG.round(k), BO_LOG.iter(k), BO_LOG.f(k), BO_LOG.L(k), BO_LOG.Lref(k), ...
            BO_LOG.w(k,1), BO_LOG.w(k,2), BO_LOG.w(k,3));
    end
    disp('==============================================================');
end

%% ===================== 11. 表2 =====================
idx_list = round(linspace(1, N, 5));
pose_list = q_actual_t(idx_list, :);

disp(' ');
disp('=============================================================================================');
disp('表2 机械臂末端位姿 (基于 URDF 齐次变换正解)');
disp('---------------------------------------------------------------------------------------------');
fprintf('序号\t px\t\t py\t\t pz\t\t psi\t theta\t phi\t\t 4x4 齐次变换矩阵 T6^0\n');
disp('---------------------------------------------------------------------------------------------');

for i = 1:size(pose_list, 1)
    q_current = pose_list(i, :);
    T = kinematics_URDF_model(q_current, p, ax);
    pos = T(1:3, 4)';

    sy = sqrt(T(1,1)^2 + T(2,1)^2);
    if sy > 1e-6
        eul = [atan2(T(2,1), T(1,1)), atan2(-T(3,1), sy), atan2(T(3,2), T(3,3))];
    else
        eul = [0, atan2(-T(3,1), sy), atan2(-T(2,3), T(2,2))];
    end

    fprintf('位姿 %d\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t\t [ %6.2f  %6.2f  %6.2f  %6.2f ]\n', ...
        i, pos(1), pos(2), pos(3), eul(3), eul(2), eul(1), T(1,:));
    fprintf('\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t [ %6.2f  %6.2f  %6.2f  %6.2f ]\n', T(2,:));
    fprintf('\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t [ %6.2f  %6.2f  %6.2f  %6.2f ]\n', T(3,:));
    fprintf('\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t [ %6.2f  %6.2f  %6.2f  %6.2f ]\n\n', T(4,:));
end

%% ===================== 12. 表3 =====================
disp('=============================================================================================');
disp('表3 逆运动学模型计算的关节变量');
disp('---------------------------------------------------------------------------------------------');
for i = 1:size(pose_list, 1)
    q_current = pose_list(i, :);
    T = kinematics_URDF_model(q_current, p, ax);
    fprintf('位姿 %d\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t\t [ %6.2f  %6.2f  %6.2f  %6.2f ]\n', ...
        i, q_current(1), q_current(2), q_current(3), q_current(4), q_current(5), q_current(6), T(1,:));
    fprintf('\t\t\t\t\t\t\t\t\t\t\t\t\t\t [ %6.2f  %6.2f  %6.2f  %6.2f ]\n', T(2,:));
    fprintf('\t\t\t\t\t\t\t\t\t\t\t\t\t\t [ %6.2f  %6.2f  %6.2f  %6.2f ]\n', T(3,:));
    fprintf('\t\t\t\t\t\t\t\t\t\t\t\t\t\t [ %6.2f  %6.2f  %6.2f  %6.2f ]\n\n', T(4,:));
end

%% ===================== 13. 图表输出 =====================
figure('Name','图表1：最优避障轨迹关节角度变化','Position',[50,550,450,350]);
plot(t, q_actual_t*180/pi, 'LineWidth',1.5);
title('最优避障轨迹关节角度'); xlabel('时间 (s)'); ylabel('角度 (°)');
legend('J1','J2','J3','J4','J5','J6','Location','best'); grid on;

figure('Name','图表1B：最优避障轨迹关节角加速度变化','Position',[520,550,450,350]);
plot(t, ddq_actual_t*180/pi, 'LineWidth',1.2);
title('最优避障轨迹关节角加速度'); xlabel('时间 (s)'); ylabel('角加速度 (°/s^2)');
legend('J1','J2','J3','J4','J5','J6','Location','best'); grid on;

figure('Name','图表2：最优避障轨迹末端坐标变化','Position',[50,100,450,350]);
plot(t, pos_actual_t, 'LineWidth',1.5);
title('最优避障轨迹末端位置'); xlabel('时间 (s)'); ylabel('位置 (m)');
legend('X','Y','Z','Location','best'); grid on;

figure('Name','图表A：轨迹三视图投影','Position',[550,50,450,850]);
subplot(3,1,1);
plot(pos_bad_t(:,1),pos_bad_t(:,2),'r--','LineWidth',1); hold on;
plot(pos_actual_t(:,1),pos_actual_t(:,2),'b-','LineWidth',2);
plot(obs.centers(:,1),obs.centers(:,2),'kX','MarkerSize',8,'LineWidth',1.5);
title('XY 平面投影'); xlabel('X (m)'); ylabel('Y (m)'); grid on; axis equal;
legend('未优化碰撞轨迹','最优避障轨迹','障碍物中心','Location','best');

subplot(3,1,2);
plot(pos_bad_t(:,2),pos_bad_t(:,3),'r--','LineWidth',1); hold on;
plot(pos_actual_t(:,2),pos_actual_t(:,3),'g-','LineWidth',2);
plot(obs.centers(:,2),obs.centers(:,3),'kX','MarkerSize',8,'LineWidth',1.5);
title('YZ 平面投影'); xlabel('Y (m)'); ylabel('Z (m)'); grid on; axis equal;

subplot(3,1,3);
plot(pos_bad_t(:,1),pos_bad_t(:,3),'r--','LineWidth',1); hold on;
plot(pos_actual_t(:,1),pos_actual_t(:,3),'m-','LineWidth',2);
plot(obs.centers(:,1),obs.centers(:,3),'kX','MarkerSize',8,'LineWidth',1.5);
title('XZ 平面投影'); xlabel('X (m)'); ylabel('Z (m)'); grid on; axis equal;

%% ===================== 14. 3D动画 =====================
h_fig3 = figure('Name','图表3：在线BO最优避障动画（IK邻域回退）','Position',[1050,200,750,700]);
ax_3d = axes('Parent',h_fig3);
plot3(ax_3d,pos_bad_t(:,1),pos_bad_t(:,2),pos_bad_t(:,3),'r--','LineWidth',1.5,'DisplayName','未优化碰撞轨迹'); hold(ax_3d,'on');
plot3(ax_3d,pos_actual_t(:,1),pos_actual_t(:,2),pos_actual_t(:,3),'b-','LineWidth',2.5,'DisplayName','在线BO最优轨迹');

[sx,sy,sz] = sphere(30);
for k = 1:obs.num
    surf(ax_3d, sx*obs.radii(k)+obs.centers(k,1), sy*obs.radii(k)+obs.centers(k,2), sz*obs.radii(k)+obs.centers(k,3), ...
        'FaceColor','r','EdgeColor','none','FaceAlpha',0.35,'DisplayName',sprintf('障碍物%d',k));
end
plot3(ax_3d,p_target_used(1),p_target_used(2),p_target_used(3),'kp','MarkerSize',12,'MarkerFaceColor','y','DisplayName','目标点');

title(ax_3d,'在线BO权重优化：3D 避障轨迹对比');
xlabel(ax_3d,'X轴(m)'); ylabel(ax_3d,'Y轴(m)'); zlabel(ax_3d,'Z轴(m)');
grid(ax_3d,'on'); axis(ax_3d,'equal'); view(ax_3d,135,30);
xlim(ax_3d,[-0.4 0.4]); ylim(ax_3d,[-0.4 0.4]); zlim(ax_3d,[0 0.6]);
legend(ax_3d,'Location','northeast');

disp('正在全速播放 3D 轨迹动画...');
anim_fps = 60;
step_skip = max(round(1/(t_sample*anim_fps)),30);
show(robot,q_actual_t(1,:),'Parent',ax_3d,'Frames','off','Visuals','off');
for i = 1:step_skip:N
    if ~isgraphics(ax_3d), break; end
    show(robot,q_actual_t(i,:),'Parent',ax_3d,'Frames','off','Visuals','off','FastUpdate',true,'PreservePlot',false);
    drawnow;
end
disp('仿真动画播放完毕！');

%% ================================ 函数区 ================================
function y = normalize01(x)
xmin = min(x); xmax = max(x);
if abs(xmax-xmin)<1e-12, y = zeros(size(x)); else, y = (x-xmin)/(xmax-xmin); end
end

function [q_sol, ok, info] = ik_solve_with_rbt(robot, q_seed, tformTarget)
persistent IK_OBJ
if isempty(IK_OBJ)
    IK_OBJ = inverseKinematics('RigidBodyTree', robot);
end
weights = [0.1 0.1 0.1 1 1 1];

seeds = [
    q_seed;
    q_seed + [ 0.25  0.15 -0.15  0.08 -0.08  0.15];
    q_seed + [-0.25 -0.15  0.15 -0.08  0.08 -0.15];
    q_seed + [-0.15 -0.35  0.25 -0.12  0.12 -0.10];
    q_seed + [ 0.70  0.10 -0.40  0.20  0.00  0.25];
    q_seed + [-0.70 -0.10  0.40 -0.20  0.00 -0.25]
];
seeds = max(-pi,min(pi,seeds));

ok = false; bestErr = inf; q_sol = q_seed; bestStatus = '';
for i = 1:size(seeds,1)
    [q_i, solInfo] = IK_OBJ('link6', tformTarget, weights, seeds(i,:));
    T_i = getTransform(robot, q_i, 'link6');
    posErr = norm(T_i(1:3,4)-tformTarget(1:3,4));

    if posErr < bestErr
        bestErr = posErr; q_sol = q_i; bestStatus = solInfo.Status;
    end
    if strcmpi(solInfo.Status,'success') && posErr < 1e-3
        ok = true; q_sol = q_i; bestErr = posErr; bestStatus = solInfo.Status; break;
    end
end
if ~ok && bestErr < 2e-3, ok = true; end
info.posErr = bestErr; info.status = bestStatus;
end

function [q_best, p_used, ok, info] = ik_fallback_nearby_search(robot, q_seed, p_target, R_target, opt)
r = opt.range; step = opt.step; gridv = -r:step:r;

cand = [];
if opt.keepSameZFirst
    for dx = gridv
        for dy = gridv
            cand = [cand; p_target + [dx,dy,0]]; %#ok<AGROW>
        end
    end
end
for dx = gridv
    for dy = gridv
        for dz = gridv
            cand = [cand; p_target + [dx,dy,dz]]; %#ok<AGROW>
        end
    end
end
cand = unique(round(cand,6),'rows');

ok = false; q_best = q_seed; p_used = p_target; bestDp = inf; bestPosErr = inf;

for i = 1:size(cand,1)
    p_try = cand(i,:);
    T_try = eye(4); T_try(1:3,1:3)=R_target; T_try(1:3,4)=p_try(:);

    [q_i, ok_i, info_i] = ik_solve_with_rbt(robot, q_seed, T_try);
    if ~ok_i, continue; end

    dp = norm(p_try - p_target); posErr = info_i.posErr;
    better = (dp < bestDp - 1e-12) || (abs(dp-bestDp)<=1e-12 && posErr < bestPosErr);

    if better
        bestDp = dp; bestPosErr = posErr; q_best = q_i; p_used = p_try; ok = true;
        if dp < 1e-12 && posErr < 1e-3, break; end
    end
end
info.bestDp = bestDp; info.bestPosErr = bestPosErr;
end

function stop = bo_output_logger(results, state)
stop = false;
if strcmp(state,'iteration')
    BO_LOG = getappdata(0,'BO_LOG');
    BO_TMP = getappdata(0,'BO_TMP');
    if isempty(BO_TMP), BO_TMP = struct('L',[],'Lref',[]); end

    i = results.NumObjectiveEvaluations;
    if i >= 1
        X = results.XTrace(i,:);
        fval = results.ObjectiveTrace(i);
        w = [10^(X.x1),10^(X.x2),10^(X.x3)];

        L = NaN; Lref = NaN;
        if length(BO_TMP.L) >= i
            L = BO_TMP.L(i);
            Lref = BO_TMP.Lref(i);
        end

        rid = getappdata(0,'BO_ROUND_ID'); if isempty(rid), rid = 0; end
        BO_LOG.round(end+1,1)=rid;
        BO_LOG.iter(end+1,1)=i;
        BO_LOG.f(end+1,1)=fval;
        BO_LOG.L(end+1,1)=L;
        BO_LOG.Lref(end+1,1)=Lref;
        BO_LOG.w(end+1,:)=w;
        setappdata(0,'BO_LOG',BO_LOG);
    end
end
end

function f = bo_outer_objective_twin(X, q_start, q_end, t_total, p, ax, obs, q_via_guess, twin)
w = [10^(X.x1),10^(X.x2),10^(X.x3)];

inner_options = optimset('Display','off','MaxIter',60,'MaxFunEvals',220,'TolX',2e-2,'TolFun',1e-3);
q_via_opt = fminsearch(@(qv) cost_function_weighted_twin(qv,q_start,q_end,t_total,p,ax,obs,w,twin), q_via_guess, inner_options);

m = evaluate_metrics_twin(q_via_opt,q_start,q_end,t_total,p,ax,obs,twin);

switch twin.online_mode
    case 'safety_balanced', alpha = [0.25,0.20,0.45,0.10];
    case 'speed_priority',  alpha = [0.20,0.15,0.55,0.10];
    case 'energy_priority', alpha = [0.20,0.15,0.35,0.30];
    otherwise,              alpha = [0.25,0.20,0.45,0.10];
end

Jacc = m.J_acc/(m.J_acc+1);
Jobs = m.J_obs/(m.J_obs+1);
Jlen = m.J_len/(m.J_len+1);
Jeng = m.J_energy/(m.J_energy+1);
Jgnd = m.J_ground/(m.J_ground+1);
Jlim = m.J_joint_limit/(m.J_joint_limit+1);

P_safe = 0;
if m.d_min_surface < twin.safety_margin
    P_safe = 2e3*(twin.safety_margin-m.d_min_surface)^2;
end

f = alpha(1)*Jacc + alpha(2)*Jobs + alpha(3)*Jlen + alpha(4)*Jeng + 0.25*Jgnd + 0.25*Jlim + P_safe;

BO_TMP = getappdata(0,'BO_TMP');
if isempty(BO_TMP), BO_TMP = struct('L',[],'Lref',[]); end
BO_TMP.L(end+1,1) = m.L_withObs;
BO_TMP.Lref(end+1,1) = m.L_ref_noObs;
setappdata(0,'BO_TMP',BO_TMP);
end

function J = cost_function_weighted_twin(q_via,q_start,q_end,t_total,p,ax,obs,w,twin)
m = evaluate_metrics_twin(q_via,q_start,q_end,t_total,p,ax,obs,twin);
J = w(1)*m.J_acc + w(2)*m.J_obs + w(3)*m.J_len + 0.05*m.J_energy + 1.0*m.J_ground + 1.0*m.J_joint_limit;
end

function m = evaluate_metrics_twin(q_via,q_start,q_end,t_total,p,ax,obs,twin)
t_eval = linspace(0,t_total,20);
[q_eval,dq_eval,ddq_eval] = trajectory_SDPOA(q_start,q_via,q_end,t_total,t_eval);

m.J_acc = sum(sum(ddq_eval.^2));
m.J_energy = sum(sum((dq_eval.^2).*twin.friction_scale_j));

q_abs = abs(q_eval);
lim_violation = max(0,q_abs - twin.joint_limit);
m.J_joint_limit = twin.joint_limit_penalty_k * sum(lim_violation(:).^2);

ee = zeros(length(t_eval),3);
J_obs = 0; J_ground = 0; d_min_surface = inf;

for i = 1:length(t_eval)
    q_i = q_eval(i,:);
    T_i = kinematics_URDF_model(q_i,p,ax);
    ee(i,:) = T_i(1:3,4)';
    ee_i = ee(i,:);

    for k = 1:obs.num
        c = obs.centers(k,:); r = obs.radii(k);
        d_center = norm(ee_i-c);
        d_surface = d_center-r;
        d_min_surface = min(d_min_surface,d_surface);

        d_safe = r + twin.safety_margin;
        margin = d_center - d_safe;
        J_obs = J_obs + log(1+exp(-18*margin));
        if d_center < d_safe
            J_obs = J_obs + 800*(d_safe-d_center)^2;
        end
    end

    jointPos = forward_all_joint_positions(q_i,p,ax);
    z_violation = max(0,twin.ground_clearance - jointPos(:,3));
    J_ground = J_ground + twin.ground_penalty_k * sum(z_violation.^2);
end

seg = diff(ee,1,1);
L_withObs = sum(sqrt(sum(seg.^2,2)));
L_ref_noObs = norm(ee(end,:)-ee(1,:));
J_len = L_withObs + 8*max(0,L_withObs - 1.2*L_ref_noObs)^2;

m.J_obs = J_obs;
m.J_len = J_len;
m.J_ground = J_ground;
m.L_withObs = L_withObs;
m.L_ref_noObs = L_ref_noObs;
m.d_min_surface = d_min_surface;
end

function jointPos = forward_all_joint_positions(q,p,ax)
T = eye(4); jointPos = zeros(6,3);
for i = 1:6
    T_trans = eye(4); T_trans(1:3,4)=p(i,:)';
    u = ax(i,:); u = u/norm(u);
    ux=u(1); uy=u(2); uz=u(3);
    th=q(i); c=cos(th); s=sin(th); v=1-c;

    R = [ux^2*v+c,      ux*uy*v-uz*s,  ux*uz*v+uy*s;
         ux*uy*v+uz*s,  uy^2*v+c,      uy*uz*v-ux*s;
         ux*uz*v-uy*s,  uy*uz*v+ux*s,  uz^2*v+c];

    T_rot = eye(4); T_rot(1:3,1:3)=R;
    T = T*(T_trans*T_rot);
    jointPos(i,:) = T(1:3,4)';
end
end

function [q_meas_t,dq_meas_t] = simulate_real_plant(q_cmd_t,dq_cmd_t,plant,dt)
N = size(q_cmd_t,1);
q_meas_t = q_cmd_t; dq_meas_t = dq_cmd_t;
for i = 1:N
    tt = (i-1)*dt;
    q_meas_t(i,:) = q_cmd_t(i,:) + plant.bias + plant.drift_rate*tt + plant.noise_std*randn(1,6);
end
dq_meas_t(2:end,:) = diff(q_meas_t,1,1)/dt;
dq_meas_t(1,:) = dq_meas_t(2,:);
d = max(0,plant.delay_steps);
if d>0 && d<N
    q_meas_t = [repmat(q_meas_t(1,:),d,1); q_meas_t(1:end-d,:)];
    dq_meas_t = [repmat(dq_meas_t(1,:),d,1); dq_meas_t(1:end-d,:)];
end
end

function err = sync_real_data(q_twin,dq_twin,q_real,dq_real)
dq = q_twin-q_real;
ddq = dq_twin-dq_real;
err.q_rmse_global = sqrt(mean(dq(:).^2));
err.dq_rmse_global = sqrt(mean(ddq(:).^2));
err.q_rmse_joint = sqrt(mean(dq.^2,1));
err.dq_rmse_joint = sqrt(mean(ddq.^2,1));
end

function twin = update_twin_params(twin,err)
g = 0.08; lambda = 0.05; ref_dq = 0.03;
e = err.dq_rmse_joint - ref_dq;
scale = max(1e-6, median(abs(e)));
delta = g*(e/scale);
twin.friction_scale_j = twin.friction_scale_j + delta - lambda*(twin.friction_scale_j-1.0);
twin.friction_scale_j = min(twin.friction_max, max(twin.friction_min, twin.friction_scale_j));
end

function T_end = kinematics_URDF_model(q,p,ax)
T_end = eye(4);
for i = 1:6
    T_trans = eye(4); T_trans(1:3,4)=p(i,:)';
    u = ax(i,:); u = u/norm(u);
    ux=u(1); uy=u(2); uz=u(3);
    th=q(i); c=cos(th); s=sin(th); v=1-c;

    R = [ux^2*v+c,      ux*uy*v-uz*s,  ux*uz*v+uy*s;
         ux*uy*v+uz*s,  uy^2*v+c,      uy*uz*v-ux*s;
         ux*uz*v-uy*s,  uy*uz*v+ux*s,  uz^2*v+c];

    T_rot = eye(4); T_rot(1:3,1:3)=R;
    T_end = T_end*(T_trans*T_rot);
end
end

function [q_ref,dq_ref,ddq_ref] = trajectory_SDPOA(q_start,q_via,q_end,t_total,t)
N = length(t);
q_ref = zeros(N,6); dq_ref = zeros(N,6); ddq_ref = zeros(N,6);

t_m = t_total/2; t_f = t_total;
M = [1,0,0,0,0,0,0;
     0,1,0,0,0,0,0;
     0,0,2,0,0,0,0;
     1,t_m,t_m^2,t_m^3,t_m^4,t_m^5,t_m^6;
     1,t_f,t_f^2,t_f^3,t_f^4,t_f^5,t_f^6;
     0,1,2*t_f,3*t_f^2,4*t_f^3,5*t_f^4,6*t_f^5;
     0,0,2,6*t_f,12*t_f^2,20*t_f^3,30*t_f^4];

for j = 1:6
    B = [q_start(j);0;0;q_via(j);q_end(j);0;0];
    C = M\B;
    q_ref(:,j) = C(1)+C(2)*t+C(3)*t.^2+C(4)*t.^3+C(5)*t.^4+C(6)*t.^5+C(7)*t.^6;
    dq_ref(:,j)= C(2)+2*C(3)*t+3*C(4)*t.^2+4*C(5)*t.^3+5*C(6)*t.^4+6*C(7)*t.^5;
    ddq_ref(:,j)=2*C(3)+6*C(4)*t+12*C(5)*t.^2+20*C(6)*t.^3+30*C(7)*t.^4;
end
end






 % %% ===================== 15. Simulink 实际物理轨迹的 FK 解算与对比 =====================
 % disp('正在进行 Simulink 实际轨迹的末端正运动学解算...');
 % 
 % % 1. 整理 Simulink 输出的数据
 % % 假设你在工作区已经提取了 6 个关节的实际输出并拼接为 actual_pos_matrix
 % % 以及对应的时间向量 sim_t (如果在 To Workspace 选了 Array，可以从 out.tout 获取)
 % % 这里我们使用你在前面图中设置的 12 个输出模块的命名逻辑：
 % sim_t = out.tout;
 % sim_theta1 = out.simout;
 % sim_theta2 = out.simout2;
 % sim_theta3 = out.simout4;
 % sim_theta4 = out.simout6;
 % sim_theta5 = out.simout8;
 % sim_theta6 = out.simout10;
 % 
 % actual_pos_matrix = [sim_theta1, sim_theta2, sim_theta3, sim_theta4, sim_theta5, sim_theta6];
 % N_sim = size(actual_pos_matrix, 1);
 % sim_end_pos = zeros(N_sim, 3);
 % 
 % % 2. 逐点进行正运动学 (Forward Kinematics) 解算
 % % 直接复用你脚本里写好的 kinematics_URDF_model 函数
 % for i = 1:N_sim
 %     q_current = actual_pos_matrix(i, :);
 %     % 计算当前时刻的 4x4 齐次变换矩阵 T_end
 %     T_current = kinematics_URDF_model(q_current, p, ax);
 %     % 提取位置向量 [px, py, pz] 并存入矩阵
 %     sim_end_pos(i, :) = T_current(1:3, 4)';
 % end
 % disp('正运动学解算完成！');
 % 
 % % 3. 绘制 3D 轨迹重合度对比图
 % figure('Name','图表4：数字孪生验证 - 理论轨迹 vs 物理仿真轨迹','Position',[150,150,800,600]);
 % 
 % % 画出之前规划好的纯数学理论轨迹 (粗蓝线)
 % plot3(pos_actual_t(:,1), pos_actual_t(:,2), pos_actual_t(:,3), 'b-', 'LineWidth', 3, 'DisplayName', '数学理论轨迹 (MATLAB)'); hold on;
 % 
 % % 画出经过电机物理约束后实际跑出来的轨迹 (红色虚线)
 % plot3(sim_end_pos(:,1), sim_end_pos(:,2), sim_end_pos(:,3), 'r--', 'LineWidth', 1.5, 'DisplayName', '物理执行轨迹 (Simulink)');
 % 
 % % 画出障碍物和目标点做参照
 % [sx, sy, sz] = sphere(30);
 % for k = 1:obs.num
 %     surf(sx*obs.radii(k)+obs.centers(k,1), sy*obs.radii(k)+obs.centers(k,2), sz*obs.radii(k)+obs.centers(k,3), ...
 %         'FaceColor','k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
 % end
 % plot3(p_target_used(1), p_target_used(2), p_target_used(3), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'DisplayName', '目标点');
 % 
 % title('机械臂末端空间位姿：控制算法与物理引擎验证');
 % xlabel('X轴 (m)'); ylabel('Y轴 (m)'); zlabel('Z轴 (m)');
 % legend('Location', 'best');
 % grid on; axis equal; view(135, 30);
 % 
 % % 4. 计算并绘制末端绝对跟踪误差 (选做)
 % % 为了计算误差，需要将 Simulink 的时间轴插值对齐到理论时间轴 t
 % sim_pos_interp = zeros(length(t), 3);
 % sim_pos_interp(:,1) = interp1(sim_t, sim_end_pos(:,1), t, 'spline');
 % sim_pos_interp(:,2) = interp1(sim_t, sim_end_pos(:,2), t, 'spline');
 % sim_pos_interp(:,3) = interp1(sim_t, sim_end_pos(:,3), t, 'spline');
 % 
 % % 计算三维空间中的欧式距离误差 (单位：毫米)
 % xyz_error_mm = sqrt(sum((pos_actual_t - sim_pos_interp).^2, 2)) * 1000;
 % 
 % figure('Name','图表5：末端空间跟踪误差','Position',[950,150,500,350]);
 % plot(t, xyz_error_mm, 'k-', 'LineWidth', 1.5);
 % title('末端笛卡尔空间综合跟踪误差');
 % xlabel('时间 (s)'); ylabel('绝对误差 (mm)');
 % grid on;

target_xyz_ts = timeseries(pos_actual_t, t);