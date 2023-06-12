% Galin Compressor Single Cylinder Simulation of Compression Stroke
%clear all;
clear torque work_em work pressure temperature time theta V Fp v x a


steel_density = 7800; %[km/m3]
vane_angle = 40; %[deg]
r = 0.041528; %[m]
V_initial = mLtom3(900);

suction_pressure = 100e3; %196e3; %[Pa]
suction_temperature = 295; %[K]
discharge_pressure = 689476;%1739.86e3; %689476;%1.6e6; %1.67e6; %[Pa], verify this (1e6 Pa == 145 psi; 689476 Pa == 100 psi; 1.2e6 == 174 psi)
polytropic_idx = 1.3;
compression_ratio = exp(log(discharge_pressure/suction_pressure) / polytropic_idx);
discharge_temperature = (discharge_pressure /suction_pressure) * (1/compression_ratio) * suction_temperature; %[K]

R = 3*r; %[m]
d = 2*r; %[m]
area = (R-r)*d; %[m^2]
lever_arm = (R+r)/2; %[m]

vane_moment_of_inertia = (2 * pi * r * steel_density * (R^4 - r^4)) * 2 * deg2rad(vane_angle) / pi;
vane_volume = (0.5*deg2rad(vane_angle)*(R^2-r^2) * 2 + pi*r^2)*d; %factor of two, because we have two vane either side of the shaft
vane_mass = vane_volume * steel_density;

V_final = V_initial / compression_ratio; %[m^3]

start_angle = V_initial / (R^2 - r^2) / d * 2; %[rad]
end_angle = V_final / (R^2 - r^2) / d * 2; %[rad]
ave_distance_travelled = lever_arm * (start_angle - end_angle); %[m]

angle_travelled = (start_angle - end_angle) / 2; %[rad] each vane travels half the angle necessary to compress from V_initial to V_final
work_done = (suction_pressure * V_initial) / (polytropic_idx - 1) * (1 - compression_ratio^(polytropic_idx-1));
torque_to_apply =  abs(work_done/2) / angle_travelled; %each vane does half the work needed to compress the gas
torque_factor = 1.00; 
torque_to_apply = torque_factor * torque_to_apply;

% simulation of just the compression stroke:
% first initialise the system, at time t=0:
idx = 1;
pressure(idx) = suction_pressure;
temperature(idx) = suction_temperature;
V(idx) = V_initial;
Fp(idx) = pressure(idx) * area; 

torque_from_gases(idx) = abs(Fp(idx)) * lever_arm;
torque_from_gases(idx) = torque_from_gases(idx); 

vane_torque_leading(idx) = -torque_to_apply + torque_from_gases(idx);
vane_torque_following(idx) = +torque_to_apply - torque_from_gases(idx);
vane_accel_leading(idx) = vane_torque_leading(idx) / vane_moment_of_inertia;
vane_accel_following(idx) = vane_torque_following(idx) / vane_moment_of_inertia;

x(idx) = start_angle;
work(idx) = 0;
vane_work_leading(idx) = 0;
vane_work_following(idx) = 0;
time(idx) = 0;

vane_pos_leading(idx) = start_angle/2;
vane_pos_following(idx) = -start_angle/2;
vane_vel_leading(idx) = 0;%152;%148;
vane_vel_following(idx) = 0;%152;%148;
bisector_pos(idx) = 0;

idx = 2;
del_t = 1e-7; %[s]
%anti-clockwise direction is positive
for t = del_t:del_t:1.0 %some arbitrary end time, long enough for model to reach V_final
    delta_speed_leading = vane_accel_leading(idx-1) * del_t;
    delta_speed_following = vane_accel_following(idx-1) * del_t;
    
    vane_vel_leading(idx) = vane_vel_leading(idx-1) + delta_speed_leading;
    vane_vel_following(idx) = vane_vel_following(idx-1) + delta_speed_following;
    
    delta_x_leading = vane_vel_leading(idx-1)*del_t + vane_accel_leading(idx-1)*del_t^2 / 2;
    delta_x_following = vane_vel_following(idx-1)*del_t + vane_accel_following(idx-1)*del_t^2 / 2;
    
    vane_pos_leading(idx) = vane_pos_leading(idx-1) + delta_x_leading;
    vane_pos_following(idx) = vane_pos_following(idx-1) + delta_x_following;
    x(idx) = vane_pos_leading(idx) - vane_pos_following(idx);
    
    time(idx) = t;
    
    V(idx) = d * (R^2 - r^2) * x(idx) * 0.5;
    
    pressure(idx) = pressure(1) * (V(1)/V(idx))^polytropic_idx;
    temperature(idx) = temperature(1)*(V(idx)/V(1))^(1-polytropic_idx);
    Fp(idx) = pressure(idx-1) * area;
    torque_from_gases(idx) = Fp(idx) * lever_arm;
    
    vane_torque_leading(idx) = -torque_to_apply + torque_from_gases(idx);
    vane_torque_following(idx) = +torque_to_apply - torque_from_gases(idx);
    vane_accel_leading(idx) = vane_torque_leading(idx) / vane_moment_of_inertia;
    vane_accel_following(idx) = vane_torque_following(idx) / vane_moment_of_inertia;

    work(idx) = torque_from_gases(idx) * delta_x_leading - torque_from_gases(idx) * delta_x_following; %this is a product of torque applied by gases, and delta_x
    vane_work_leading(idx) = -torque_to_apply * delta_x_leading;
    vane_work_following(idx) = +torque_to_apply * delta_x_following;
    
    bisector_pos(idx) = (vane_pos_leading(idx)+vane_pos_following(idx))/2;
    if (V(idx) <= V_initial/compression_ratio)
        fprintf("ended loop as max compression achieved\n");
        break;
    end
    if (bisector_pos(idx) >= pi/2) 
        fprintf("ended loop as reached 90 deg rev\n");
        break;
    end
    
    idx = idx+1;
end

fprintf("galin engine - compression stroke stats: \n");
fprintf("sanity check values:\n");
fprintf("compression ratio []: %f\n",compression_ratio);
fprintf("pressure difference, start to end [kPa]: %f -> %f [%f -> %f (psi)]\n", pressure(1)/1e3, pressure(end)/1e3, pressure(1)*0.000145038, pressure(end)*0.000145038);
fprintf("volume difference, start to end [mL]: %f -> %f\n", V(1)*1e6, V(end)*1e6);
fprintf("temperature difference, start to end [K]: %f -> %f\n", temperature(1), temperature(end));


fprintf("work done on gases during compression [J]: %f\n", sum(work));
fprintf("work done by leading vane [J]: %f\n", sum(vane_work_leading));
fprintf("work done by following vane [J]: %f\n", sum(vane_work_following));
fprintf("stroke time [ms]: %f\n",time(end)*1000);
fprintf("average torque leading [Nm]: %f\n",mean(vane_torque_leading));
fprintf("average torque following [Nm]: %f\n",mean(vane_torque_following));
fprintf("max magnitude of velocity reached [rad/s]: %f\n",abs(min(vane_vel_leading)));
fprintf("discharge pressure [MPa]: %f ( %f psi) \n",pressure(end)/1e6,pressure(end)*0.000145038);
fprintf("discharge temperature [K]: %f\n",temperature(end));
fprintf("velocity of rotation [RPM]: %f\n", 1/(4*time(end))*60);
fprintf("wanted bisector speed estimate [rad/s]: %f\n", (pi/2)/time(end));
fprintf("torque requested from electrical machines [Nm]: %f\n",torque_to_apply);
fprintf("electrical machine power req. (each) [kW]: %f\n", torque_to_apply*(pi/2/time(end))/1000)
fprintf("CFM: %f\n",V_initial*1000*8*(1/(4*time(end)))*2.1188799727597);
fprintf("compression ratio: %f\n", V_initial/V(idx));
fprintf("leading vane location [deg]: %f\n", vane_pos_leading(end)*360/pi/2);
fprintf("rotation of bisector [deg]: %f\n", ((vane_pos_leading(end) + vane_pos_following(end)) / 2)*360/2/pi);
fprintf("CFM/HP: %f\n", (V_initial*1000*8*(1/(4*time(end)))*2.1188799727597)/(torque_to_apply*(pi/2/time(end))/1000*2));