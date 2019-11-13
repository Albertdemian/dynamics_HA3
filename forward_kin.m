clear all
close all
clc
%% 

% links = struct('a1' , 'd1')
links.a2 = 4;
links.d1 = 10; 

% joints = struct('t1' , 't2','d3')
joints.t1 = 20;
joints.t2 = 10;
joints.d3 = 1;
%%

t1 = deg2rad(joints.t1);
t2 = deg2rad(joints.t2);
d3 = deg2rad(joints.d3);

a2 = links.a2;
d1 = links.d1; 

%% DH-Parameter
syms th1 th2 d3 
alp1 = deg2rad(90);
alp2 = deg2rad(-90);
alp3 = deg2rad(0);


transmission.t01 = [th1 d1 0 alp1];
transmission.t12 = [th2 0 0 alp2];
transmission.t23 = [0  d3+a2 0 alp3];
transmission.notation = ["theta" "alpha"]

%% Rotation matrices 

[rotat trans transform03] = rotation(transmission);

%% Forward Kinematics 
th1 = t1;
th2 = t2;
d3 = joints.d3
eff_pos = eval(trans.o03)
position.Px = eff_pos(1);
position.Py = eff_pos(2);
position.Pz = eff_pos(3);
position.position = eff_pos'

%% Inverse Kinematics 
r = sqrt(position.Px^2 + position.Py^2);
s = position.Pz - d1;

inv.theta1 =180+ rad2deg(atan2(position.Py,position.Px)); 
inv.theta2 =90- rad2deg(atan2(s,r));
inv.d3  = sqrt(r^2+s^2) - a2; 
inv.inverse = [inv.theta1 inv.theta2 inv.d3]



%% Jacobian Matrix
%Geometrical approach 

v01 = cross(([0 0 1]'),trans.o03); 
v02 = cross((rotat.r01*[0 0 1]'),(trans.o03-trans.o01));
v03 = rotat.r02*[0 0 1]';
w01 = rotat.r01*[0 0 1]';
w02 = rotat.r02*[0 0 1]';
w03 = [0 0 0]';
J_geo = eval( [v01 v02 v03;w01 w02 w03])

% Analytical approach

O = transform03(1:3,4)
Z = transform03(1:3,3);
OZ = [O ; Z]
V1 = diff(OZ,'th1')
V2 = diff(OZ,'th2')
V3 = diff(OZ,'d3')
J_an = [V1 V2 V3]
J_an_re = eval(J_an)


%% Singularities 

J = J_an(1:3,1:3);
deter = det(J)==0
sing = solve(deter,th1,th2,d3)

%% path  
index = 1;
for i = 0:0.10:10
th1 = sin(i);
th2 = cos(2*i);
d3 = sin(3*i);
transmission.t01 = [th1 d1 0 alp1];
transmission.t12 = [th2 0 0 alp2];
transmission.t23 = [0  d3+a2 0 alp3];
[ro, position, trans]= rotation(transmission);

pos_list(index,:)= position.o03';
index = index+1;
th1_list(index)=th1;
th2_list(index)=th2;
d3_list(index)=d3;
end

figure()
subplot(2,2,1)
plot3(pos_list(:,1),pos_list(:,2),pos_list(:,3))
grid on
title('3d path')
xlabel('x position')
ylabel('y position')
zlabel('z position')
legend('path in world frame')

subplot(2,2,2)
plot(pos_list(:,1),pos_list(:,2))
grid on
title('X-Y projection of path')
xlabel('x position')
ylabel('y position')

subplot(2,2,3)
plot(pos_list(:,2),pos_list(:,3))
grid on
title('Y-Z projection of path')
xlabel('Y position')
ylabel('Z position')

subplot(2,2,4)
plot(pos_list(:,1),pos_list(:,3))
grid on
title('X-Z projection of path')
xlabel('X position')
ylabel('Z position')

figure('Name','trajectory')
plot(th1_list,'r')
hold on
plot(th2_list,'g')
hold on 
plot(d3_list,'b')
legend('theta 1','theta 2','d3')
title('trajectory')

%% Trajectory
index =1;

for i = 0:0.1:10
    position.Px = 2*a2*sin(i);
    position.Py = 2*a2*sin(2*i);
    position.Pz = d1 * sin(3*i);
    position.position = eff_pos';

    r = sqrt(position.Px^2 + position.Py^2);
    s = position.Pz - d1;

    inv.theta1 =180+ rad2deg(atan2(position.Py,position.Px)); 
    inv.theta2 =90- rad2deg(atan2(s,r));
    inv.d3  = sqrt(r^2+s^2) - a2; 
    inv.inverse = [inv.theta1 inv.theta2 inv.d3]
       
    path_x(index) = position.Px;
    path_y(index) = position.Py;
    path_z(index) = position.Pz;
    trajectory(index,:)=inv.inverse;
   
    index = index+1;
end
dvar(1,1:3)=0;
for i = 2:length(trajectory)
dvar(i,:) = (trajectory(i,:)-trajectory(i-1,:))/0.1;
end
t = 0:0.1:10
figure('Name','trajectory plan')
subplot(2,2,1)
plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3))
grid on
title('3d path')
xlabel('theta 1')
ylabel('theta 2')
zlabel('d3')
legend('trajectory in joint space')

subplot(2,2,2)
plot(t,dvar(:,1))
grid on
title('th1 velocity')
xlabel('time')
ylabel('theta 1')

subplot(2,2,3)
plot(t,dvar(:,2))
grid on
title('theta 2 velocity')
xlabel('time')
ylabel('theta 2')

subplot(2,2,4)
plot(t,dvar(:,3))
grid on
title('d3 velocity')
xlabel('time')
ylabel('d 3')

%% solving for speed
th1_vel = 5;
th2_vel = 5;
d3_vel  = 4;

th1_acc = 1; 
th2_acc = 1;
d3_acc  =1; 

tau = max([th1_vel th2_vel d3_vel] ./ [th1_acc th2_acc d3_acc])

T = max((trajectory(2,:)-trajectory(1,:))./[th1_vel th2_vel d3_vel])

if (th1_vel > th1_acc*tau) || (th2_vel > th2_acc*tau) || (d3_vel > d3_acc*tau)
    th1_acc = th1_vel/tau;
    th2_acc = th2_vel/tau;
    d3_acc  = d3_vel/tau;
end


accel_t = 0:0.1:tau;
v_accel = [th1_acc*accel_t' th2_acc*accel_t' d3_acc*accel_t'];
const_t = tau:0.1 :T 
v_const = v_accel(end,:).*ones(length(const_t),1)
dec = T:0.1:T+tau
dec_t = tau:-0.1:0
v_dec = [th1_acc*dec_t' th2_acc*dec_t' d3_acc*dec_t']

tajectory_v = [v_accel;v_const;v_dec]
time_traj = [accel_t'; const_t'; dec']
figure()
subplot(3,1,1)
plot(time_traj, tajectory_v(:,1))
title('theta 1 trajectory')
xlabel('t')
ylabel('theta1_dot')

subplot(3,1,2)
plot(time_traj, tajectory_v(:,2))
title('theta 2 trajectory')
xlabel('t')
ylabel('theta2_dot')

subplot(3,1,3)
plot(time_traj, tajectory_v(:,3))
title('d3 trajectory')
xlabel('t')
ylabel('d3_dot')