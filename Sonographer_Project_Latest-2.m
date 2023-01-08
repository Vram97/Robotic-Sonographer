%% Sonographer Project 
% This script analyzes a 4-DOF robot designed to apply an ultrasound probe
% at three locations while applying 3 lbs of force to the body.
%% Setup
%clearing all existing variable values and closing all windows
clear all
close all
counter = 1;
%Defining values
P1 = [-175; 125; 0]; %mm
P2 = [175; 125; 0]; %mm
P3 = [0; 150; -200]; %mm
F1 = [0.961; -0.961; 0]; %kg
F2 = [-0.961; -0.961; 0]; %kg
F3 = [0; -1.36; 0]; %kg
% applying initial conditions (all joint variables set to 0)
q = [0, 0, 0, 0];
% applying initial conditions (new_approximation being set to values not
% equal to joint variables)
new_approximation = [1, 1, 1, 1];
%defining length values (meters)
L1 = 75; %length of link 1
L2 = 350; %length of link 2
L3 = 400; %length of link 3
L4 = 150; %length of link 4
L_body = 400; %distance between robot space frame and origin of body in the x direction
%defining mass values
m_motor = 2.4/2.2; %kg
m_butterfly = 0.309; %kg
%calculating the masses of the links using the volume from Autodesk Fusion
%360 and the density of aluminum, and adding the attached components
d_aluminum = 0.098; %lb/in^3
vol_1 = 8.282; %in^3
vol_2 = 16.003; %in^3
vol_3 = 18.138; %in^3
vol_4 = 6.945; %in^3
m1 = d_aluminum*vol_1/2.2 + m_motor; % weight of link plus weight of motor, kg
m2 = d_aluminum*vol_2/2.2 + m_motor; % weight of link plus weight of motor, kg
m3 = d_aluminum*vol_3/2.2 + m_motor; % weight of link plus weight of motor, kg
m4 = d_aluminum*vol_4/2.2 + m_motor + m_butterfly; % weight of link plus weight of motor and ultrasound, kg

%% Forward and Inverse Kinematics
% This section determines the forward and inverse kinematics of the robot
% to setting the loop to continuously run the newton-raphson method for
% inverse kinematics until the joint variables are equal to the new 
% approximation values and the variables converge
for k = 0:1:2
if k == 0
    q = [0, pi/2, -pi/2, -pi/4];
end
if k == 1
    q = [0, pi/4, -pi/4, -pi*.75];
end
if k == 2
    q = [pi/4, pi/2, -pi/2, -pi/2];
end
i = 0;
while ((q(1) ~= new_approximation(1)) || (q(2) ~= new_approximation(2))) && (q(3) ~= new_approximation(3))
%after the first run, set q to the calculated "new approximation" joint
%value for each iteration of the newton-raphson method
if i == 1
    q = new_approximation';
end
%set i = 1 after the initial if statement to set q = new_approximation
%after the first run of this while loop
i = 1;

%M matrix for the home position of the robot
M = [0, 0,  1, L2 + L3 + L4;
     0, -1, 0, 0;
     1, 0, 0, L1;
     0, 0, 0,  1];

%w and v screw axes and velocity vectors for each joint, followed by
%calculating the skewomega of each joint, the rotation matrix of each
%joint, the distance vector for each joint, then the exponential function
%of each joint
w1 = [0; 0; 1];
v1 = [0; 0; 0];
skewomega1 = [0, -w1(3), w1(2); w1(3), 0, -w1(1); -w1(2), w1(1), 0];
R1 = eye(3,3) + sin(q(1))*skewomega1 + (1-cos(q(1)))*skewomega1^2;
d1 =  (eye(3,3)*q(1) + (1-cos(q(1)))*skewomega1 + (q(1)-sin(q(1)))*skewomega1^2)*v1;
eSTheta1 = [R1, d1;
    0, 0, 0, 1];
w2 = [0; -1; 0];
v2 = [L1; 0; 0];
skewomega2 = [0, -w2(3), w2(2); w2(3), 0, -w2(1); -w2(2), w2(1), 0];
R2 = eye(3,3) + sin(q(2))*skewomega2 + (1-cos(q(2)))*skewomega2^2;
d2 = (eye(3,3)*q(2) + (1-cos(q(2)))*skewomega2 + (q(2)-sin(q(2)))*skewomega2^2)*v2;
eSTheta2 = [R2, d2;
    0, 0, 0, 1];
w3 = [0; -1; 0];
v3 = [L1; 0; -L2];
skewomega3 = [0, -w3(3), w3(2); w3(3), 0, -w3(1); -w3(2), w3(1), 0];
R3 = eye(3,3) + sin(q(3))*skewomega3 + (1-cos(q(3)))*skewomega3^2;
d3 =  (eye(3,3)*q(3) + (1-cos(q(3)))*skewomega3 + (q(3)-sin(q(3)))*skewomega3^2)*v3;
eSTheta3 = [R3, d3;
    0, 0, 0, 1];
w4 = [0; -1; 0];
v4 = [L1; 0; -(L2+L3)];
skewomega4 = [0, -w4(3), w4(2); w4(3), 0, -w4(1); -w4(2), w4(1), 0];
R4 = eye(3,3) + sin(q(4))*skewomega4 + (1-cos(q(4)))*skewomega4^2;
d4 =  (eye(3,3)*q(4) + (1-cos(q(4)))*skewomega4 + (q(4)-sin(q(4)))*skewomega4^2)*v4;
eSTheta4 = [R4, d4;
    0, 0, 0, 1];

%calculating the product of the exponentials to get the Transformation
%Matrix from the base joint to the end effector
T = eSTheta1*eSTheta2*eSTheta3*eSTheta4*M;
if counter == 1
T_first = T;
end
counter = 2;
%% 
%defining other transformation matrices to plot robot while converging to
%the desired location
M1 = eye(4, 4);
M2 = [1, 0, 0, 0;
      0, 0, -1, 0;
      0, 1, 0, L1;
      0, 0, 0, 1];
M3 = [1, 0, 0, L2;
      0, 0, -1, 0;
      0, 1, 0, L1;
      0, 0, 0, 1];
M4 = [1, 0, 0, L2+L3;
      0, 0, -1, 0;
      0, 1, 0, L1;
      0, 0, 0, 1];
M5 = M;
T01 = M1;
T02 = eSTheta1*M2;
T03 = eSTheta1*eSTheta2*M3;
T04 = eSTheta1*eSTheta2*eSTheta3*M4;
T05 = T;

%plotting robot while converging to the desired location
figure (1)
axis on
hold on
view([1,-1,1])
plot3([0, T01(1,4)],...
      [0, T01(2,4)],...
      [0, T01(3,4)], '-k', 'LineWidth', 3)

plot3([T01(1,4), T02(1,4)],...
      [T01(2,4), T02(2,4)],...
      [T01(3,4), T02(3,4)], '-o', 'LineWidth', 3)

plot3([T02(1,4), T03(1,4)],...
      [T02(2,4), T03(2,4)],...
      [T02(3,4), T03(3,4)], '-k', 'LineWidth', 3)

plot3([T03(1,4), T04(1,4)],...
      [T03(2,4), T04(2,4)],...
      [T03(3,4), T04(3,4)], '-o', 'LineWidth', 3)

plot3([T04(1,4), T05(1,4)],...
      [T04(2,4), T05(2,4)],...
      [T04(3,4), T05(3,4)], '-k', 'LineWidth', 3)

%Calculating the Adjoints of each individual step needed to calculate the
%Jacobian vectors.
[Adt_1] = adjoint(eSTheta1);
[Adt_2] = adjoint(eSTheta1*eSTheta2);
[Adt_3] = adjoint(eSTheta1*eSTheta2*eSTheta3);
[Adt_4] = adjoint(eSTheta1*eSTheta2*eSTheta3*eSTheta4);

%Defining the S matrix to multiply columns of it by the Adjoints calculated
%previously. This matrix is the combination of the w and v vectors defined
%within the previous section of code.
S = [w1 w2 w3 w4;
     v1 v2 v3 v4];

%Manipulating the S matrix and multiplying sections of it (where necessary)
%to the previously calculated Adjoints to calculate the space Jacobian.
V1 = S(:,1);
V2 = Adt_1*S(:,2);
V3 = Adt_2*S(:,3);
V4 = Adt_3*S(:,4);

%Calculating the space Jacobian
J_s = [V1, V2, V3, V4];

%Converting the space Jacobian to the Body Jacobian.
J_b = adjoint(inv(T))*J_s;

%Calculating the Analytical Jacobian by taking the linear velocity portion
%of the Body Jacobian and multiplying the rotational portion of the
%transformation matrix from the space frame to the body frame by it.
J_v = J_b(4:6, :);
R_sb = T(1:3, 1:3);
J_analytical = R_sb*J_v;


%Defining the homogeneous transformation matrix T_sd and calculating T_bd,
%the transformation matrix in the body frame for the desired pose
T_sd1 = [-sqrt(2)/2 0 sqrt(2)/2 L_body+P1(1); 0 -1 0 P1(3); -sqrt(2)/2 0 -sqrt(2)/2 P1(2); 0 0 0 1]; %First pose
T_sd2 = [-sqrt(2)/2 0 -sqrt(2)/2 L_body+P2(1); 0 -1 0 P2(3); sqrt(2)/2 0 -sqrt(2)/2 P2(2); 0 0 0 1]; %second pose
T_sd3 = [1 0 0 L_body+P3(1); 0 -1 0 -P3(3); 0 0 -1 P3(2); 0 0 0 1]; %Third pose
T_bd1 = inv(T)*T_sd1;
T_bd2 = inv(T)*T_sd2;
T_bd3 = inv(T)*T_sd3;

%Separating the Rotation Matrix and position vector from the transformation
%matrix T_bd, and doing it for each run (position 1, 2, and 3)
if k == 0
    R_bd = T_bd1(1:3, 1:3);
    p_bd = T_bd1(1:3, 4);
end
if k == 1
    R_bd = T_bd2(1:3, 1:3);
    p_bd = T_bd2(1:3, 4);
end
if k == 2
    R_bd = T_bd3(1:3, 1:3);
    p_bd = T_bd3(1:3, 4);
end

%Calculating theta, the skew symmetric matrix of omega, the inverse G
%function of theta, and the linear velocity vector v_bd to be able to
%create the vector S
theta = acos(0.5*(trace(R_bd)-1));
skew_w = (1/(2*sin(theta)))*(R_bd - transpose(R_bd));
inv_G_theta = (1/theta)*eye(3,3) - 0.5*skew_w + ((1/theta) - 0.5*cot(theta/2))*skew_w^2;
v_bd = inv_G_theta*p_bd;
S_bd = [skew_w(3,2); skew_w(1,3); skew_w(2,1); v_bd];

%Using the vector S as well as the calculated theta to calculate the twist
%in the body frame to the desired pose from the transform T_bd
V_b = S_bd*theta;

%Approximating the iterative change in joint values to move the end 
%effector to the desired/specified location.
delta_theta = pinv(J_b)*V_b;

%Adding the calculated change to the current joint values to get the new
%approximated joint values for the first iteration of the Newton Raphson
%Method.
new_approximation = q' + delta_theta;

end
%% Plotting 
% Plotting the final robot pose after the newton-raphson method has made the robot converge to the desired location
if k == 0
    figure (2)
    axis on
    hold on
    view([1,-1,1])
    plot3([0, T01(1,4)],...
        [0, T01(2,4)],...
        [0, T01(3,4)], '-k', 'LineWidth', 3)

    plot3([T01(1,4), T02(1,4)],...
        [T01(2,4), T02(2,4)],...
        [T01(3,4), T02(3,4)], '-o', 'LineWidth', 3)

    plot3([T02(1,4), T03(1,4)],...
        [T02(2,4), T03(2,4)],...
        [T02(3,4), T03(3,4)], '-k', 'LineWidth', 3)

    plot3([T03(1,4), T04(1,4)],...
        [T03(2,4), T04(2,4)],...
        [T03(3,4), T04(3,4)], '-o', 'LineWidth', 3)

    plot3([T04(1,4), T05(1,4)],...
        [T04(2,4), T05(2,4)],...
        [T04(3,4), T05(3,4)], '-k', 'LineWidth', 3)
    figure (5)
    axis on
    hold on
    view([1,-1,1])
    plot3([0, T01(1,4)],...
        [0, T01(2,4)],...
        [0, T01(3,4)], '-k', 'LineWidth', 3)

    plot3([T01(1,4), T02(1,4)],...
        [T01(2,4), T02(2,4)],...
        [T01(3,4), T02(3,4)], '-o', 'LineWidth', 3)

    plot3([T02(1,4), T03(1,4)],...
        [T02(2,4), T03(2,4)],...
        [T02(3,4), T03(3,4)], '-k', 'LineWidth', 3)

    plot3([T03(1,4), T04(1,4)],...
        [T03(2,4), T04(2,4)],...
        [T03(3,4), T04(3,4)], '-o', 'LineWidth', 3)

    plot3([T04(1,4), T05(1,4)],...
        [T04(2,4), T05(2,4)],...
        [T04(3,4), T05(3,4)], '-k', 'LineWidth', 3)
end
if k == 1
    figure (3)
    axis on
    hold on
    view([1,-1,1])
    plot3([0, T01(1,4)],...
        [0, T01(2,4)],...
        [0, T01(3,4)], '-k', 'LineWidth', 3)

    plot3([T01(1,4), T02(1,4)],...
        [T01(2,4), T02(2,4)],...
        [T01(3,4), T02(3,4)], '-o', 'LineWidth', 3)

    plot3([T02(1,4), T03(1,4)],...
        [T02(2,4), T03(2,4)],...
        [T02(3,4), T03(3,4)], '-k', 'LineWidth', 3)

    plot3([T03(1,4), T04(1,4)],...
        [T03(2,4), T04(2,4)],...
        [T03(3,4), T04(3,4)], '-o', 'LineWidth', 3)

    plot3([T04(1,4), T05(1,4)],...
        [T04(2,4), T05(2,4)],...
        [T04(3,4), T05(3,4)], '-k', 'LineWidth', 3)
    figure (5)
    axis on
    hold on
    view([1,-1,1])
    plot3([0, T01(1,4)],...
        [0, T01(2,4)],...
        [0, T01(3,4)], '-k', 'LineWidth', 3)

    plot3([T01(1,4), T02(1,4)],...
        [T01(2,4), T02(2,4)],...
        [T01(3,4), T02(3,4)], '-o', 'LineWidth', 3)

    plot3([T02(1,4), T03(1,4)],...
        [T02(2,4), T03(2,4)],...
        [T02(3,4), T03(3,4)], '-k', 'LineWidth', 3)

    plot3([T03(1,4), T04(1,4)],...
        [T03(2,4), T04(2,4)],...
        [T03(3,4), T04(3,4)], '-o', 'LineWidth', 3)

    plot3([T04(1,4), T05(1,4)],...
        [T04(2,4), T05(2,4)],...
        [T04(3,4), T05(3,4)], '-k', 'LineWidth', 3)
end
if k == 2
    figure (4)
    axis on
    hold on
    view([1,-1,1])
    plot3([0, T01(1,4)],...
        [0, T01(2,4)],...
        [0, T01(3,4)], '-k', 'LineWidth', 3)

    plot3([T01(1,4), T02(1,4)],...
        [T01(2,4), T02(2,4)],...
        [T01(3,4), T02(3,4)], '-o', 'LineWidth', 3)

    plot3([T02(1,4), T03(1,4)],...
        [T02(2,4), T03(2,4)],...
        [T02(3,4), T03(3,4)], '-k', 'LineWidth', 3)

    plot3([T03(1,4), T04(1,4)],...
        [T03(2,4), T04(2,4)],...
        [T03(3,4), T04(3,4)], '-o', 'LineWidth', 3)

    plot3([T04(1,4), T05(1,4)],...
        [T04(2,4), T05(2,4)],...
        [T04(3,4), T05(3,4)], '-k', 'LineWidth', 3)
    figure (5)
    axis on
    hold on
    view([1,-1,1])
    plot3([0, T01(1,4)],...
        [0, T01(2,4)],...
        [0, T01(3,4)], '-k', 'LineWidth', 3)

    plot3([T01(1,4), T02(1,4)],...
        [T01(2,4), T02(2,4)],...
        [T01(3,4), T02(3,4)], '-o', 'LineWidth', 3)

    plot3([T02(1,4), T03(1,4)],...
        [T02(2,4), T03(2,4)],...
        [T02(3,4), T03(3,4)], '-k', 'LineWidth', 3)

    plot3([T03(1,4), T04(1,4)],...
        [T03(2,4), T04(2,4)],...
        [T03(3,4), T04(3,4)], '-o', 'LineWidth', 3)

    plot3([T04(1,4), T05(1,4)],...
        [T04(2,4), T05(2,4)],...
        [T04(3,4), T05(3,4)], '-k', 'LineWidth', 3)
end
% outputting the transformation matrix and joint values when the
% newton-raphson method converges to the solution
    T_solved = T
    q_solved = q
      
%% Joint Torque Calculation
%defining theta, Vi0, and Vi_dot0
theta = (q);
Vi0 = zeros(6,1);
Vi_dot0 = [0,0,0,0,0,9.81]';

%defining the constants
theta_1_dot = 0; 
theta_2_dot = 0;
theta_3_dot = 0;
theta_4_dot = 0;
theta_1_double = 0; 
theta_2_double = 0;
theta_3_double = 0;
theta_4_double = 0;

%defining the transformation matrices between each frame when the robot is
%in the home configuration
M01 = eye(4,4);

M12 = [1, 0, 0, 0;
       0, 0, -1, 0;
       0, 1, 0, L1;
       0, 0, 0, 1];

M23 = [1, 0, 0, L2;
       0, 1, 0, 0;
       0, 0, 1, 0;
       0, 0, 0, 1];

M34 = [1, 0, 0, L3;
       0, 1, 0, 0;
       0, 0, 1, 0;
       0, 0, 0, 1];

M45 = [0, 0, 1, L4;
       1, 0, 0, 0;
       0, 1, 0, 0;
       0, 0, 0, 1];

%Defining the M and A matrices for the Newton Euler method
M(:,:,1) = M1;
M(:,:,2) = M2;
M(:,:,3) = M3;
M(:,:,4) = M4;
M(:,:,5) = M5;

numjoints = 4;
A = zeros(6, numjoints);
 for l = 1:numjoints
     A(:,l) = adjoint2(S(:,l),inv(M(:,:,l)));
 end

%Calculating the transformation matrices between each joint when the robot
%is in the desired pose
T_01 = M1 * twist2ht(S(:,1), theta(1));
T10 = inv(T_01);

T_12 = M12 * twist2ht(S(:,2), theta(2));
T21 = inv(T_12);

T_23 = M23 * twist2ht(S(:,3), theta(3));
T32 = inv(T_23);

T_34 = M34 * twist2ht(S(:,4), theta(4));
T43 = inv(T_34);

T54 = [0, 0, 1, L4;
       1, 0, 0, 0;
       0, 1, 0, 0;
       0, 0, 0, 1];

%calculating the velocities and accelerations of each joint
V1 = adjoint2(Vi0,T10) + A(:,1)*theta_1_dot;
V2 = adjoint2(V1,T21) + A(:,2)*theta_2_dot;
V3 = adjoint2(V2,T32) + A(:,3)*theta_3_dot;
V4 = adjoint2(V3,T43) + A(:,4)*theta_4_dot;
Vd1 = adjoint2(Vi_dot0,T10) + ad(V1)*A(:,1)*theta_1_dot + A(:,1)*theta_1_double;
Vd2 = adjoint2(Vd1,T21) + ad(V2)*A(:,2)*theta_2_dot + A(:,2)*theta_2_double;
Vd3 = adjoint2(Vd2,T32) + ad(V3)*A(:,3)*theta_3_dot + A(:,3)*theta_3_double;
Vd4 = adjoint2(Vd3,T43) + ad(V4)*A(:,4)*theta_4_dot + A(:,4)*theta_4_double;

%defining the spatial inertia matrices
G1 = [zeros(3,3), zeros(3,3);
      zeros(3,3), m1*eye(3,3)];
G2 = [zeros(3,3), zeros(3,3);
      zeros(3,3), m2*eye(3,3)];
G3 = [zeros(3,3), zeros(3,3);
      zeros(3,3), m3*eye(3,3)];
G4 = [zeros(3,3), zeros(3,3);
      zeros(3,3), m4*eye(3,3)];

%Defining Ftip for each case and calculating the forces and torques
%required by each joint to apply the desired force at the end effector
if k == 0
Ftip = [F1(1) F1(2) F1(3) 0 0 0]';
end
if k == 1
Ftip = [F2(1) F2(2) F2(3) 0 0 0]';
end
if k == 2
Ftip = [F3(1) F3(2) F3(3) 0 0 0]';
end
Fi_4 = adjoint(T54)'*Ftip + G4*Vd4 - ad(V4)'*G4*V4;
Fi_3 = adjoint(T43)'*Fi_4 + G3*Vd3 - ad(V3)'*G3*V3;
Fi_2 = adjoint(T32)'*Fi_3 + G2*Vd2 - ad(V2)'*G2*V2;
Fi_1 = adjoint(T21)'*Fi_2 + G1*Vd1 - ad(V1)'*G1*V1;
tau4 = Fi_4'*A(:,4)/1000 %N*m
tau3 = Fi_3'*A(:,3)/1000 %N*m
tau2 = Fi_2'*A(:,2)/1000 %N*m
tau1 = Fi_1'*A(:,1)/1000 %N*m
end

%% Functions Used
%The functions used in this code are defined below:

function [Adt] = adjoint(T)
%get rotation matrix and vector d from T
R = T(1:3, 1:3);
d = T(1:3, 4);
%calculate skew symmetric matrix for d
skew_d=[0 -d(3) d(2); d(3) 0 -d(1); -d(2) d(1) 0];
%calculate V's representation in respect to new frame T
Adt=[R zeros(3,3);skew_d*R R];
end

function Vtrans = adjoint2(V,T)
%get rotation matrix and vector d from T
R = [T(1) T(5) T(9);
    T(2) T(6) T(10);
    T(3) T(7) T(11)];
d = [T(13) T(14) T(15)];
%calculate skew symmetric matrix for d
skew_d=[0 -d(3) d(2); d(3) 0 -d(1); -d(2) d(1) 0];
%calculate V's representation in respect to new frame T
Vtrans=[R zeros(3,3);skew_d*R R]*V;
end

function adV = ad(V)
w = V(1:3);
v = V(4:6);
skew_w = [0, -w(3), w(2); w(3), 0, -w(1); -w(2), w(1), 0];
skew_v = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
adV = [skew_w zeros(3,3);skew_v skew_w];
end

function T = twist2ht(S,theta)
omega = [S(1) S(2) S(3)];
v = transpose([S(4) S(5) S(6)]);
skewomega = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
R = eye(3,3) + sin(theta)*skewomega + (1-cos(theta))*skewomega^2;
position = (eye(3,3)*theta + (1-cos(theta))*skewomega + (theta - sin(theta))*skewomega^2)*v;
T = [R, position;
    0, 0, 0, 1];
end