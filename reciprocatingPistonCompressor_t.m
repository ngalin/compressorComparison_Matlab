% Reciprocating Piston Single Cylinder Simulation of Compression Stroke
% - simulating the compression stroke in a reciprocating piston setup with
% a crank shaft is complicated by the fact that we don't know what torque
% to apply to the shaft to get the compression to occur. We guess at the
% torque value by calculating the work done to compress the gases, and then
% dividing that work by the distance the piston has to travel from Vmax to
% Vmin (i.e. 2 x crank_radius)
clear all;

% constants:
gas_constant = 8.3145; %[J/mol K]


% simulation parameters:
steel_density = 7800; %[km/m3]
vane_angle = deg2rad(40); %[deg]
r = 0.041528; %[m]
%Vc = mLtom3(100);
polytropic_idx = 1.3;

suction_pressure = 100e3; %196e3; %[Pa]
suction_temperature = 295; %[K]
discharge_pressure = 689476; %1739.86e3; %689476;%1.6e6; %1.67e6; %[Pa], verify this (1e6 Pa == 145 psi; 689476 Pa == 100 psi; 1.2e6 == 174 psi)
compression_ratio = exp(log(discharge_pressure/suction_pressure) / polytropic_idx);

R = 3*r; 
d = 2*r; 
vane_area = (R-r)*d; %[m^2]
lever_arm = (R+r)/2; %[m]
vane_moment_of_inertia = (2 * pi * r * steel_density * (R^4 - r^4)) * 2 * deg2rad(vane_angle) / pi;
vane_volume = (0.5*deg2rad(vane_angle)*(R^2-r^2) * 2 + pi*r^2)*d; %factor of two, because we have two vane either side of the shaft
vane_mass = vane_volume * steel_density;

V_initial = mLtom3(900); %[m^3]
V_final = V_initial / compression_ratio;% - V_initial; %[m^3]
Vc = V_final;

piston_mass = 2 * vane_mass;
piston_mass = piston_mass / 2; %NOTICE undoing the doubling of the mass
bore = sqrt((2 * vane_area) * 4 / pi); %doubling the area because in a RP engine the gases only act on the piston head, instead of on the areas of two vanes either side of combustion chamber
bore = bore / 2; %NOTICE undoing the doubling of the bore area
crank_radius = 2 * (V_initial - V_final) / (pi*bore^2);
stroke = 2 * crank_radius;
ratio_connecting_rod_length_2_crank_radius = 4; %heywood, p43 
connecting_rod_length = ratio_connecting_rod_length_2_crank_radius * crank_radius; %[m]
n = ratio_connecting_rod_length_2_crank_radius;
piston_area = pi * (bore/2)^2;

work_done = (suction_pressure * V_initial) / (polytropic_idx - 1) * (1 - compression_ratio^(polytropic_idx-1));
force_to_apply = work_done / (2*crank_radius);

% simulation:
% first initialise the system, at time t=0:
idx = 1;
pressure(idx) = suction_pressure;
temperature(idx) = suction_temperature;

theta(idx) = 180;
lever_arm(idx) = crank_radius * (sind(theta(idx)) + sind(2*theta(idx)) / (2*sqrt((connecting_rod_length/crank_radius)^2 - sind(theta(idx))^2)));
Fp(idx) = pressure(idx) * piston_area;
Fp(idx) = Fp(idx) + force_to_apply;
a(idx) = Fp(idx) / piston_mass;
torque(idx) = Fp(idx)*lever_arm(idx);
x(idx) = 2*crank_radius;
s(idx) = crank_radius*cosd(theta(idx)) + sqrt(connecting_rod_length^2 - crank_radius^2*sind(theta(idx))^2);
V(idx) = Vc + pi*bore^2/4*(connecting_rod_length + crank_radius - s(idx));
v(idx) = 0;
total_work(idx) = 0;
work_em(idx) = 0;
time(idx) = 0;

idx = 2;
del_t = 1e-7; %[s]
for t = del_t:del_t:10 %some arbitrary end time, long enough for model to reach condition when theta > 180 deg.    
    delta_speed = a(idx-1)*del_t;
    v(idx) = v(idx-1) + delta_speed;
    delta_x = v(idx-1)*del_t + a(idx-1)*del_t^2 / 2;
    x(idx) = x(idx-1) + delta_x;
    time(idx) = t;
    
    cos_theta = (x(idx)^2 + 2*crank_radius^2 - 2*x(idx)*crank_radius - 2*x(idx)*n*crank_radius + 2*n*crank_radius^2) / (2*crank_radius^2 + 2*n*crank_radius^2 - 2*x(idx)*crank_radius);
    theta(idx) = abs(acosd(cos_theta)); 
    
    %knowing the crank angle allows us to calculate cylinder volume:
    s(idx) = crank_radius*cosd(theta(idx)) + sqrt(connecting_rod_length^2 - crank_radius^2*sind(theta(idx))^2);
    V(idx) = Vc + pi*bore^2/4*(connecting_rod_length + crank_radius - s(idx));
    
    pressure(idx) = pressure(1) * (V(1)/V(idx))^polytropic_idx;
    temperature(idx) = temperature(1)*(V(idx)/V(1))^(1-polytropic_idx);
    Fp(idx) = pressure(idx-1) * piston_area;
    Fp(idx) = Fp(idx) + force_to_apply;
    a(idx) = Fp(idx) / piston_mass;

    lever_arm(idx) = crank_radius * (sind(theta(idx)) + sind(2*theta(idx)) / (2*sqrt((connecting_rod_length/crank_radius)^2 - sind(theta(idx))^2)));
    torque(idx) = Fp(idx) * lever_arm(idx);
    
    total_work(idx) = Fp(idx) * delta_x;
    work_em(idx) = force_to_apply * delta_x;
    
    if (theta(idx) <= 0) 
        fprintf("ended because theta got to > 180 deg\n");
        fprintf("at the end the ratio of max(V) to min(V) is: %f\n", max(V)/min(V));
        break;
    end
    
    if (V_initial/V(idx) >= (compression_ratio-2.5310e-08))
        fprintf("ended because compression ratio achieved: %f\n", max(V)/min(V));
        fprintf("at the end theta is %f\n",theta(end));
        break;
    end
    
    idx = idx+1;
end

rp_torque = torque;
rp_work = work_em;
rp_pressure = pressure;
rp_temperature = temperature;
rp_time = time;
rp_theta = theta;
rp_volume = V;
rp_force = Fp;
rp_vel = v;
rp_dist = x;
rp_accel = a;
rp_lever_arm = lever_arm;
    
%remember that this is one stroke out of 4, the other 3 don't have any
%torque or work generation.
fprintf("reciprocating piston - compression stroke stats: \n");
fprintf("sanity check values:\n");
fprintf("compression ratio []: %f\n", compression_ratio);
fprintf("work done on gases during compression [J]: %f\n", sum(total_work));
fprintf("work done by electrical machine [J]: %f\n", sum(work_em));
fprintf("pressure difference, start to end [kPa]: %f -> %f [%f -> %f (psi)]\n", pressure(1)/1e3, pressure(end)/1e3, pressure(1)*0.000145038, pressure(end)*0.000145038);
fprintf("volume difference, start to end [mL]: %f -> %f\n", V(1)*1e6, V(end)*1e6);
fprintf("temperature difference, start to end [K]: %f -> %f\n", temperature(1), temperature(end));

fprintf("output values:\n");

fprintf("stroke time [ms]: %f\n",time(end)*1000);
fprintf("based on the stroke time (time to rotate 180 deg), the revolutions per second that the electrical machine needs to make is [ms]: %f\n", 1/(time(end)*2));
fprintf("average torque [Nm]: %f\n",mean(torque));
fprintf("average speed of piston [rad/s]: %f\n", mean(v));
fprintf("torque requested from electrical machine as 0.8 of max [Nm]: %f\n", max(abs(force_to_apply*lever_arm))*0.8);
fprintf("based on torque(0.8 of max) and stroke time, power req. of electrical machine [kW]: %f\n", max(abs(torque))*0.8*(1/(time(end)*2)*2*pi)/1000);
fprintf("based on the torque load on EM machine * 0.8 we calculate the power of the EM [kW]: %f\n", max(abs(force_to_apply*lever_arm))*0.8 * (1/(time(end)*2)*2*pi) / 1000);
fprintf("CFM: %f\n",V_initial*1000*1*(1/(2*time(end)))*2.1188799727597);
fprintf("CFM/HP: %f\n", (V_initial*1000*1*(1/(2*time(end)))*2.1188799727597)/(max(torque)*0.8*(1/(time(end)*2)*2*pi)/1000));
