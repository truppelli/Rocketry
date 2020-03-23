function integrate_inu_data;
% This function is designed to test the integration of numerical data 
% saved from the intertial navigation unit gyroscope and accelerometer 
% sensors.
%
% The function starts by defining a series of processing parameters that
% are used to load the data file, extract the important information,
% perform calculations based upon this data, and output data plots.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Processing Parameters                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the filename of the inertial sensor CSV data to read.  This 
% string should contain the full pathname to the file.

data_filename_read= '\Users\Phillip Truppelli\Desktop\LOGGER01.CSV';
%'~/Desktop/LOGGER01.CSV';

% This is a Boolean value stating whether to plot the raw unprocessed
% acceleration and rotation data (this should be set as true at first to
% find the time range over which to perform the trajectory calculations).
% If this value is set to true, then only the raw data is plotted.  If this
% is set to false, then the data is processed and the full data is plotted.
plot_raw_data=false;

% This is a Boolean value stating whether to use the built in Euler angles
% (if True) or to use the integrated Euler angles (if False)
use_sensor_euler_angles=false;

% This is a Boolean value stating whether to smooth the data or not
smoothing_on=true;
% This is the type of smoothing filter to apply to the data
smoothing_filter_type='sgolay';
% This is the span (ie the number of points 'averaged' over)
smoothing_filter_span=10;

% This is the time range over which to perform the calculations (in
% seconds)
t_min=1000.482;
t_max=1150.17;

% This specifies the interpolation algorithm method to use for
% interpolating the sensor data
interp_method='cubic';

% This specifies the finite difference method used to estimate the time
% derivative of the rotational velocity
finite_diff_method='4th_order_richardson_noise_min';

% This is the defined initial condition of the position, ie this is
% defining the origin
r_x_0=0;
r_y_0=0;
r_z_0=0;
% This defines the initial conditions of the linear velocity, ie this is
% the velocity of the sensor at the start of the numerical integration
v_x_0=0;
v_y_0=0;
v_z_0=0;
% This defines the initial conditions of the angular position, ie this is
% the angle of the sensor at the start of the numerical integration
theta_x_0=0;
theta_y_0=0;
theta_z_0=0;

% These are the constant offsets for the angular velocity (ie if there is a
% positive bias error - this should be the negative of that value)
omega_x_constant=0;
omega_y_constant=0;
omega_z_constant=0;
% These are the scale factors for the angular velocity (ie if the measured
% angular velocity is scaled by some factor from the actual angular
% velocity, this is that scale factor)
omega_x_scale=1;
omega_y_scale=1;
omega_z_scale=1;

% These are the constant offsets for the acceleration (ie if there is a
% positive bias error - this should be the negative of that value)
a_x_constant=0;
a_y_constant=0;
a_z_constant=0;
% These are the scale factors for the acceleration (ie if the measured
% acceleration is scaled by some factor from the actual acceleration,
% this is that scale factor)
a_x_scale=1;
a_y_scale=1;
a_z_scale=1;

% This defines the figure window size
window_size=[500,350];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading the Data File                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This reads in the CSV inertial sensor data and saves it into the
% 'inertial_data' structure
inertial_data=read_csv_inu_data(data_filename_read);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the Raw Data (if Specified in the Parameters Section)          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% This sets the units of the root object (the screen) to pixels
set(0,'units','pixels');
% This extracts the screen resolution
screen_size=get(0,'screensize');

% This closes any open figure windows
close all;

% If the user defines that the raw data should be plotted, this plots it
if plot_raw_data;
    
    % This specifies the time range over which to extract the data (and
    % sets it to the full time range)
    t_index_min=1;
    t_index_max=length(inertial_data.data.omega_x);
    
    % This extracts the time vector
    t=inertial_data.data.t(t_index_min:t_index_max)/1e3;
    
    % This extracts the angular velocity data vectors
    omega_x=omega_x_scale*inertial_data.data.omega_x(t_index_min:t_index_max)+omega_x_constant;
    omega_y=omega_y_scale*inertial_data.data.omega_y(t_index_min:t_index_max)+omega_y_constant;
    omega_z=omega_z_scale*inertial_data.data.omega_z(t_index_min:t_index_max)+omega_z_constant;
    
    % This extracts the linear acceleration vectors
    a_lin_x=a_x_scale*inertial_data.data.a_lin_x(t_index_min:t_index_max)+a_x_constant;
    a_lin_y=a_y_scale*inertial_data.data.a_lin_y(t_index_min:t_index_max)+a_y_constant;
    a_lin_z=a_z_scale*inertial_data.data.a_lin_z(t_index_min:t_index_max)+a_z_constant;
    
    % This plots the rotational velocity data
    h1=figure(1);
    set(h1,'OuterPosition',[1,screen_size(4)-window_size(2),window_size(1)+1,window_size(2)+1]);
    plot(t,omega_x,'-RED','Linewidth',2);
    hold on;
    plot(t,omega_y,'-BLUE','Linewidth',2);
    plot(t,omega_z,'-MAGENTA','Linewidth',2);
    hold off;
    grid on;
    xlabel('Time [s]','FontSize',18);
    ylabel('\Omega [1/s]','FontSize',18);
    legend('\Omega_x','\Omega_y','\Omega_z');
    title('Rotational Velocity versus Time','FontSize',18);
    set(gca,'FontSize',18);
    
    % This plots the linear acceleration data
    h2=figure(2);
    set(h2,'OuterPosition',[window_size(1)+2,screen_size(4)-window_size(2),window_size(1)+1,window_size(2)+1]);
    plot(t,a_lin_x,'-RED','Linewidth',2);
    hold on;
    plot(t,a_lin_y,'-BLUE','Linewidth',2);
    plot(t,a_lin_z,'-MAGENTA','Linewidth',2);
    hold off;
    grid on;
    xlabel('Time [s]','FontSize',18);
    ylabel('a [m/s^2]','FontSize',18);
    legend('a_x','a_y','a_z');
    title('Acceleration versus Time','FontSize',18);
    set(gca,'FontSize',18);
    
    % This returns and quits the script
    return;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessig the Inertial Sensor Data                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This calculates the indices which as closely represent the time interval
% listed
t_index_min=interp1(inertial_data.data.t/1e3,1:length(inertial_data.data.t),t_min,'nearest');
t_index_max=interp1(inertial_data.data.t/1e3,1:length(inertial_data.data.t),t_max,'nearest');

% This extracts the time vector
t=inertial_data.data.t(t_index_min:t_index_max)/1e3;

% This extracts the angular velocity data vectors
omega_x=omega_x_scale*inertial_data.data.omega_x(t_index_min:t_index_max)+omega_x_constant;
omega_y=omega_y_scale*inertial_data.data.omega_y(t_index_min:t_index_max)+omega_y_constant;
omega_z=omega_z_scale*inertial_data.data.omega_z(t_index_min:t_index_max)+omega_z_constant;

% If specified by the user, this performs a smoothing operation on the
% rotational velocity data
if (smoothing_on);
    % This smooths the rotational velocity data in each dimension
    omega_x_smooth=smooth(omega_x,smoothing_filter_span,smoothing_filter_type);
    omega_y_smooth=smooth(omega_y,smoothing_filter_span,smoothing_filter_type);
    omega_z_smooth=smooth(omega_z,smoothing_filter_span,smoothing_filter_type);
else;
    % This copies the rotational velocity data into the "smooth" variables
    % for easier coding
    omega_x_smooth=omega_x;
    omega_y_smooth=omega_y;
    omega_z_smooth=omega_z;
end;

if use_sensor_euler_angles;
    
    % This extracts the Euler angles from the sensor data
    theta_x=inertial_data.data.theta_x(t_index_min:t_index_max);
    theta_y=inertial_data.data.theta_y(t_index_min:t_index_max);
    theta_z=inertial_data.data.theta_z(t_index_min:t_index_max);
    
    % This converts the Euler angles from degrees to radians
    theta_x=theta_x*pi/180;
    theta_y=theta_y*pi/180;
    theta_z=theta_z*pi/180;
    
    % This unwraps the angle data
    theta_x=unwrap(theta_x,2*pi);
    theta_y=unwrap(theta_y,2*pi);
    theta_z=unwrap(theta_z,2*pi);
    
else;
    
    % This performs the numerical integration on the rotational velocity data
    % to calculate the rotational position versus time
    theta_x=cumulative_uneven_simpsons_integration(t,omega_x_smooth);
    theta_y=cumulative_uneven_simpsons_integration(t,omega_y_smooth);
    theta_z=cumulative_uneven_simpsons_integration(t,omega_z_smooth);

end;

% This sets the angles equal to the initial conditions at the start of the
% integration time
theta_x=theta_x-theta_x(1)+theta_x_0;
theta_y=theta_y-theta_y(1)+theta_y_0;
theta_z=theta_z-theta_z(1)+theta_z_0;

% This converts the Euler angles to quaternions
rotation_quaternion=angle2quat(theta_x,theta_y,theta_z);

% This extracts the linear acceleration vectors
a_lin_x=a_x_scale*inertial_data.data.a_lin_x(t_index_min:t_index_max)+a_x_constant;
a_lin_y=a_y_scale*inertial_data.data.a_lin_y(t_index_min:t_index_max)+a_y_constant;
a_lin_z=a_z_scale*inertial_data.data.a_lin_z(t_index_min:t_index_max)+a_z_constant;

% If specified by the user, this performs a smoothing operation on the
% linear acceleration data
if (smoothing_on);
    % This smooths the rotational velocity data in each dimension
    a_lin_x_smooth=smooth(a_lin_x,smoothing_filter_span,smoothing_filter_type);
    a_lin_y_smooth=smooth(a_lin_y,smoothing_filter_span,smoothing_filter_type);
    a_lin_z_smooth=smooth(a_lin_z,smoothing_filter_span,smoothing_filter_type);
else;
    % This copies the rotational velocity data into the "smooth" variables
    % for easier coding
    a_lin_x_smooth=a_lin_x;
    a_lin_y_smooth=a_lin_y;
    a_lin_z_smooth=a_lin_z;
end;

% This creates a composite array of all the linear acceleration vectors
% throughout the extracted time period
a_lin_smooth=[zeros(size(a_lin_x)),a_lin_x_smooth,a_lin_y_smooth,a_lin_z_smooth];

% This rotates the linear acceleration vectors into the world coordinate
% system (at t=t_min) by multiplying the acceleration vector by the
% quaternion rotation
a_lin_world=quatmultiply(quatmultiply(rotation_quaternion,a_lin_smooth),quatconj(rotation_quaternion));

% This extracts the individual vector components out from the world
% coordinate linear acceleration vector
a_lin_x_world=a_lin_world(:,2);
a_lin_y_world=a_lin_world(:,3);
a_lin_z_world=a_lin_world(:,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically Integrating Acceleration ODE                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This creates data interpolant functions of the linear acceleration, the
% rotational velocity, and the rotational acceleration
data_interpolant=create_data_interpolants(t,a_lin_x_world,a_lin_y_world,a_lin_z_world,...
    omega_x_smooth,omega_y_smooth,omega_z_smooth,...
    interp_method,finite_diff_method);

% This defines the initial conditions of the differential equations
% governing the equations of motion
X0=[v_x_0;r_x_0;v_y_0;r_y_0;v_z_0;r_z_0];

% This numerically integrates the equations of motion over the specified
% time interval
[t_ode,X]=ode45(@(T,X)rotational_acceleration_ode(T,X,data_interpolant),[min(t),max(t)],X0);

% This extracts the numerical integrated position and velocity vectors
v_x_ode=X(:,1);
r_x_ode=X(:,2);
v_y_ode=X(:,3);
r_y_ode=X(:,4);
v_z_ode=X(:,5);
r_z_ode=X(:,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Calculated Trajectory Data                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This plots the rotational velocity data
h1=figure(1);
set(h1,'OuterPosition',[1,screen_size(4)-window_size(2),window_size(1)+1,window_size(2)+1]);
plot(t,omega_x_smooth,'-RED','Linewidth',2);
hold on;
plot(t,omega_y_smooth,'-BLUE','Linewidth',2);
plot(t,omega_z_smooth,'-MAGENTA','Linewidth',2);
hold off;
grid on;
xlabel('Time [s]','FontSize',18);
ylabel('\Omega [1/s]','FontSize',18);
legend('\Omega_x','\Omega_y','\Omega_z');
title('Rotational Velocity versus Time','FontSize',18);
set(gca,'FontSize',18);

% This plots the rotational position data
h2=figure(2);
set(h2,'OuterPosition',[window_size(1)+2,screen_size(4)-window_size(2),window_size(1)+1,window_size(2)+1]);
plot(t,mod(theta_x,2*pi),'-RED','Linewidth',2);
hold on;
plot(t,mod(theta_y,2*pi),'-BLUE','Linewidth',2);
plot(t,mod(theta_z,2*pi),'-MAGENTA','Linewidth',2);
hold off;
grid on;
xlabel('Time [s]','FontSize',18);
ylabel('\theta [rad]','FontSize',18);
legend('\theta_x','\theta_y','\theta_z');
title('Rotational Position versus Time','FontSize',18);
set(gca,'FontSize',18);

% This plots the linear acceleration data
h3=figure(3);
set(h3,'OuterPosition',[1,screen_size(4)-2*window_size(2),window_size(1)+1,window_size(2)+1]);
plot(t,a_lin_x_world,'-RED','Linewidth',2);
hold on;
plot(t,a_lin_y_world,'-BLUE','Linewidth',2);
plot(t,a_lin_z_world,'-MAGENTA','Linewidth',2);
hold off;
grid on;
xlabel('Time [s]','FontSize',18);
ylabel('a [m/s^2]','FontSize',18);
legend('a_x','a_y','a_z');
title('Acceleration (in World Coordinates) versus Time','FontSize',18);
set(gca,'FontSize',18);

% This plots the linear velocity data
h4=figure(4);
set(h4,'OuterPosition',[window_size(1)+2,screen_size(4)-2*window_size(2),window_size(1)+1,window_size(2)+1]);
plot(t_ode,v_x_ode,'-RED','Linewidth',2);
hold on;
plot(t_ode,v_y_ode,'-BLUE','Linewidth',2);
plot(t_ode,v_z_ode,'-MAGENTA','Linewidth',2);
hold off;
grid on;
xlabel('Time [s]','FontSize',18);
ylabel('v [m/s]','FontSize',18);
legend('v_x','v_y','v_z');
title('Velocity (in World Coordinates) versus Time','FontSize',18);
set(gca,'FontSize',18);

% This plots the linear position data
h5=figure(5);
set(h5,'OuterPosition',[2*window_size(1)+2,screen_size(4)-2*window_size(2),window_size(1)+1,window_size(2)+1]);
plot(t_ode,r_x_ode,'-RED','Linewidth',2);
hold on;
plot(t_ode,r_y_ode,'-BLUE','Linewidth',2);
plot(t_ode,r_z_ode,'-MAGENTA','Linewidth',2);
hold off;
grid on;
xlabel('Time [s]','FontSize',18);
ylabel('r [m]','FontSize',18);
legend('r_x','r_y','r_z');
title('Position (in World Coordinates) versus Time','FontSize',18);
set(gca,'FontSize',18);

% This plots the 3D trajectory data
h6=figure(6);
set(h6,'OuterPosition',[1,screen_size(4)-3*window_size(2),window_size(1)+1,window_size(2)+1]);
plot3(r_x_ode,r_y_ode,r_z_ode,'-BLACK','Linewidth',2);
axis equal;
grid on;
xlabel('x [m]','FontSize',18);
ylabel('y [m]','FontSize',18);
zlabel('z [m]','FontSize',18);
title('3D Trejectory (in World Coordinates)','FontSize',18);
set(gca,'FontSize',18);



function inertial_data=read_csv_inu_data(data_filename_read);
% This function reads in the inertial navigation unit comma seperated value
% (CSV) data file and saves the data into the structure 'interial_data'.

% This imports the data from the SD card file
imported_data=importdata(data_filename_read,',',4);
% This reads in the UNIX time data
unix_time_cell=textscan(imported_data.textdata{1,1},'%*s %*s %u');
% This extracts the start time at which the script began ranning since 
% midnight on January 1, 1970 in seconds
inertial_data.time.unix=unix_time_cell{1,1};
% This reads in the standard time string
standard_time_cell=textscan(imported_data.textdata{2,1},'%*s %*s %*["]%d%*[/]%d%*[/]%d %d%*[:]%d%*[:]%d%*["]',1);
% This extracts the standard time
inertial_data.time.year=standard_time_cell{1,1};
inertial_data.time.month=standard_time_cell{1,2};
inertial_data.time.day=standard_time_cell{1,3};
inertial_data.time.hour=standard_time_cell{1,4};
inertial_data.time.minute=standard_time_cell{1,5};
inertial_data.time.second=standard_time_cell{1,6};
% This reads in the data unit strings
data_units_cell=textscan(imported_data.textdata{3,1},'%s');
% This is the number of fields in the data
data_field_number=length(data_units_cell{1});
% This initializes a cell array of the units
inertial_data.units=cell(1,data_field_number);
% This iterates through the units saving them to memory
for unit_index=1:data_field_number;
    % This extracts the current index string
    current_unit_string=data_units_cell{1}{unit_index};
    % If this isn't the last column (which won't contain a comma), then the
    % comma is removed
    if (unit_index<data_field_number);
        % This removes the comma at the end of the unit (this could probably be
        % done with the 'textscan' call, but regexp are hard . . . )
        current_unit_string=current_unit_string(1:end-1);
    end;
    % This adds the current unit to the data unit cell array
    inertial_data.units{1,unit_index}=current_unit_string;
end;
% This initializes a cell array of the variables
inertial_data.variables=cell(1,data_field_number);
% This iterates through the variables saving them to memory
for variable_index=1:data_field_number;
    % This extracts the current index string
    current_variable_string=imported_data.colheaders{variable_index};
    % If this isn't the last column (which won't contain a comma), then the
    % comma is removed
    if (variable_index>1);
        % This removes the comma at the end of the variable (this could probably be
        % done with the 'textscan' call, but regexp are hard . . . )
        current_variable_string=current_variable_string(2:end);
    end;
    % This adds the current unit to the data unit cell array
    inertial_data.variables{1,variable_index}=current_variable_string;
    % This extracts the data for the current variable
    eval(['inertial_data.data.',current_variable_string,'=imported_data.data(:,',num2str(variable_index),');']);
end;



function dXdt=rotational_acceleration_ode(t_interp,X,data_interpolant);
% This function calculates the values of the differential equations
% governing the position of the inertial sensor unit.  In particular this
% function returns the values of the rotational differential equations.

% This extracts the position and linear velocity information from the 'X' 
% variable for calculating the equations of motion
vx=X(1);
rx=X(2);
vy=X(3);
ry=X(4);
vz=X(5);
rz=X(6);

% This calculates the interpolated values of the acceleration data at the
% current time instant
a_x_interp=ppval(data_interpolant.a_x,t_interp);
a_y_interp=ppval(data_interpolant.a_y,t_interp);
a_z_interp=ppval(data_interpolant.a_z,t_interp);

% This calculates the interpolated values of the angular velocity data at
% the current time instant
omega_x_interp=ppval(data_interpolant.omega_x,t_interp);
omega_y_interp=ppval(data_interpolant.omega_y,t_interp);
omega_z_interp=ppval(data_interpolant.omega_z,t_interp);

% This calculates the interpolated values of the angular acceleration data 
% at the current time instant
d_omega_x_dt_interp=ppval(data_interpolant.d_omega_x_dt,t_interp);
d_omega_y_dt_interp=ppval(data_interpolant.d_omega_y_dt,t_interp);
d_omega_z_dt_interp=ppval(data_interpolant.d_omega_z_dt,t_interp);

% This calculates the values of the derivatives of the equations of motion
% governing rotational acceleration
dXdt=rotational_acceleration_equations(rx,ry,rz,...
    vx,vy,vz,...
    a_x_interp,a_y_interp,a_z_interp,...
    omega_x_interp,omega_y_interp,omega_z_interp,...
    d_omega_x_dt_interp,d_omega_y_dt_interp,d_omega_z_dt_interp);



function dXdt=rotational_acceleration_equations(rx,ry,rz,vx,vy,vz,ax,ay,az,wx,wy,wz,dwxdt,dwydt,dwzdt);
% This function returns the equations of motion describing rotational
% acceleration.  In particular the function takes the input scalar 
% arguments:
%
%   rx, ry, rz              Linear position in each dimension
%   vx, vy, vz              Linear velocity in each dimension
%   ax, ay, az              Linear acceleration in each dimension
%   wx, wy, wz              Rotational velocity in each dimension
%   dwxdt, dwydt, dwzdt     Rotational acceleration in each dimension
%
% and outputs a 6x1 vector 'dXdt' that contains the differential equations
% describing the rotational acceleration.  In particular the components of
% 'dXdt' are given by:
%
%   dXdt(1)     dvxdt       Linear acceleration in x-dimension
%   dXdt(2)     drxdt       Linear velocity in x-dimension
%   dXdt(3)     dvydt       Linear acceleration in y-dimension
%   dXdt(4)     drydt       Linear velocity in y-dimension
%   dXdt(5)     dvzdt       Linear acceleration in z-dimension
%   dXdt(6)     drzdt       Linear velocity in z-dimension
%
% The output from these equations are then used as the output for a
% numerical ODE solver to determine the position versus time and velocity
% versus time of the sensor.

% This is the equation of motion describing the linear acceleration in the
% x-dimension for a rotating reference frame
dvxdt=ax+2*(wy*vz-wz*vy)+wx*wy*ry-(wy^2)*rx+wx*wz*rz-(wz^2)*rx-dwydt*rz+dwzdt*ry;
% This is the equation of motion describing the linear velocity in the
% x-dimension for a rotating reference frame
drxdt=vx;
% This is the equation of motion describing the linear acceleration in the
% y-dimension for a rotating reference frame
dvydt=ay+2*(wz*vx-wx*vz)-(wx^2)*ry+wx*wy*rx+wy*wz*rz-(wz^2)*ry+dwxdt*rz-dwzdt*rx;
% This is the equation of motion describing the linear velocity in the
% y-dimension for a rotating reference frame
drydt=vy;
% This is the equation of motion describing the linear acceleration in the
% z-dimension for a rotating reference frame
dvzdt=az+2*(wx*vy-wy*vx)-(wx^2)*rz-(wy^2)*rz+wx*wz*rx+wy*wz*ry-dwxdt*ry+dwydt*rx;
% This is the equation of motion describing the linear velocity in the
% z-dimension for a rotating reference frame
drzdt=vz;

% This initializes the output vector
dXdt=zeros(6,1);

% This writes the equations of motion derivative values to the output
% vector
dXdt(1)=dvxdt;
dXdt(2)=drxdt;
dXdt(3)=dvydt;
dXdt(4)=drydt;
dXdt(5)=dvzdt;
dXdt(6)=drzdt;



function data_interpolant=create_data_interpolants(t,a_x,a_y,a_z,omega_x,omega_y,omega_z,interp_method,finite_diff_method);
% This function creates a data structure containing interpolant functions
% for evaluating the sensor data at interpolated time steps.

% This creates an interpolant function for the acceleration data
a_x_interpolant=interp1(t,a_x,interp_method,'pp');
a_y_interpolant=interp1(t,a_y,interp_method,'pp');
a_z_interpolant=interp1(t,a_z,interp_method,'pp');

% This creates an interpolant function for the rotational velocity data
omega_x_interpolant=interp1(t,omega_x,interp_method,'pp');
omega_y_interpolant=interp1(t,omega_y,interp_method,'pp');
omega_z_interpolant=interp1(t,omega_z,interp_method,'pp');

% This calculates the average time interval between measurements, ie the
% optimal spacing to interpolate the functions over in some sense
dt_mean=mean(diff(t));

% This creates a vector of times to inerpolate the rotational velocity
% data over so that it has even spacing and finite difference methods may
% be used to estimate the rotation acceleration
t_interp_vector=linspace(min(t),max(t),round((max(t)-min(t))/dt_mean)+1);

% This interpolates the rotational velocity data onto the constantly spaced
% time intervals
omega_x_interp=interp1(t,omega_x,t_interp_vector,interp_method);
omega_y_interp=interp1(t,omega_y,t_interp_vector,interp_method);
omega_z_interp=interp1(t,omega_z,t_interp_vector,interp_method);

% This calculates the finite difference estimates of the time derivatives
% of the rotational velocity, ie an estimate of the rotational
% acceleration.  The strange formatting of the function is due to a bug in
% 'scalar_field_gradient' that doesn't let it calculate gradients on 1D
% functions - so the the code pretends that the data is 2D.
[d_omega_x_dt,~]=scalar_field_gradient([omega_x_interp',omega_x_interp'],finite_diff_method);
[d_omega_y_dt,~]=scalar_field_gradient([omega_y_interp',omega_x_interp'],finite_diff_method);
[d_omega_z_dt,~]=scalar_field_gradient([omega_z_interp',omega_x_interp'],finite_diff_method);

% Due to the strange form of the last function call to calculate the time
% derivative of the data, the time derivative vector is duplicated.  This
% removes the duplicate data and permutes the array to have consistent
% formatting.
d_omega_x_dt=(d_omega_x_dt(:,1)/dt_mean)';
d_omega_y_dt=(d_omega_y_dt(:,1)/dt_mean)';
d_omega_z_dt=(d_omega_z_dt(:,1)/dt_mean)';

% This creates an interpolant function for the rotational acceleration data
d_omega_x_dt_interpolant=interp1(t,d_omega_x_dt,interp_method,'pp');
d_omega_y_dt_interpolant=interp1(t,d_omega_y_dt,interp_method,'pp');
d_omega_z_dt_interpolant=interp1(t,d_omega_z_dt,interp_method,'pp');

% This creates a data structure to contain the interpolant functions
data_interpolant=struct();
% This adds the linear acceleration interpolant functions to the data
% structure
data_interpolant.a_x=a_x_interpolant;
data_interpolant.a_y=a_y_interpolant;
data_interpolant.a_z=a_z_interpolant;
% This adds rotational velocity interpolant functions to the data structure
data_interpolant.omega_x=omega_x_interpolant;
data_interpolant.omega_y=omega_y_interpolant;
data_interpolant.omega_z=omega_z_interpolant;
% This adds the rotational acceleration interpolant functions to the data
% structure
data_interpolant.d_omega_x_dt=d_omega_x_dt_interpolant;
data_interpolant.d_omega_y_dt=d_omega_y_dt_interpolant;
data_interpolant.d_omega_z_dt=d_omega_z_dt_interpolant;



function y_integrate=cumulative_uneven_simpsons_integration(x,y);
% This function calculates a cumulative integration of the function given 
% by the vectors 'x' and 'y' which must be of the same length.  The 
% function is assumed to be of the form
%	y(i) = f(x(i))
% where the spacing of the 'x' values need not be constant.  The function 
% uses a form of Simpson's rule to calculate a numerical approximation to 
% the integral of the function.  This method fits a quadratic function to 
% each set of three points and calculates the exact definite integral of 
% the quadratic function.  For all middle points there are two quadratic 
% fits to each pair of points, so this takes an average over both fits - 
% this roughly halves the error.

% This initializes the cumulative numerical integral
y_integrate=zeros(size(x));

% This initializes the x variable arrays to all zeros
x1=zeros(size(x));
x2=zeros(size(x));
x3=zeros(size(x));
% This sets the x variable arrays to offsets of the x coordinates
x1(3:end)=x(1:end-2);
x2(2:end)=x(1:end-1);
x3(1:end)=x(1:end);

% This initializes the y variable arrays to all zeros
y1=zeros(size(x));
y2=zeros(size(x));
y3=zeros(size(x));
% This sets the y variable arrays to offsets of the y coordinates
y1(3:end)=y(1:end-2);
y2(2:end)=y(1:end-1);
y3(1:end)=y(1:end);

% This calculates the definite integral from x1 to x2
y_quad_12=((x1-x2).*(3*(x3.^2).*(y1+y2)-2*x2.*x3.*(2*y1+y2)+(x2.^2).*(y1-y3)+(x1.^2).*(y2-y3)+2*x1.*(-x3.*(y1+2*y2)+x2.*(y1+y2+y3))))./(6*(x1-x3).*(-x2+x3));
% This calculates the definite integral from x2 to x3
y_quad_23=((x2-x3).*((x3.^2).*(y1-y2)+(x2.^2).*(y1-y3)-3*(x1.^2).*(y2+y3)+2*x1.*x3.*(2*y2+y3)-2*x2.*(x3.*(y1+y2+y3)-x1.*(y2+2*y3))))./(6*(x1-x2).*(x1-x3));

% This sets the second element of the cumulative integral equal to the
% first half of the first Simpson's integral
y_integrate(2)=y_quad_12(3);
% This sets the middle elements of the cumulative integral equal to the
% average of the first and second half of the Simpson's integrals from
% subsequent calculations
y_integrate(3:end-1)=(y_quad_12(4:end)+y_quad_23(3:end-1))/2;
% This sets the last element of the cumulative integral equal to the second
% half of the last Simpson's integral
y_integrate(end)=y_quad_23(end);

% This calculates the cumulative sum of all of the numerical integrations
y_integrate=cumsum(y_integrate);



function q = angle2quat( r1, r2, r3, varargin )
%  ANGLE2QUAT Convert rotation angles to quaternion.
%   Q = ANGLE2QUAT( R1, R2, R3 ) calculates the quaternion, Q, for given,
%   R1, R2, R3.   R1 is an M array of first rotation angles.  R2 is an M
%   array of second rotation angles.  R3 is an M array of third rotation
%   angles.  Q returns an M-by-4 matrix containing M quaternions. Q has its
%   scalar number as the first column.  Rotation angles are input in radians.    
%
%   Q = ANGLE2QUAT( R1, R2, R3, S ) calculates the quaternion, Q, for a
%   given set of rotation angles, R1, R2, R3, and a specified rotation
%   sequence, S.  
%
%   The default rotation sequence is 'ZYX' where the order of rotation
%   angles for the default rotation are R1 = Z Axis Rotation, R2 = Y Axis
%   Rotation, and R3 = X Axis Rotation. 
%
%   All rotation sequences, S, are supported: 'ZYX', 'ZYZ', 'ZXY', 'ZXZ',
%   'YXZ', 'YXY', 'YZX', 'YZY', 'XYZ', 'XYX', 'XZY', and 'XZX'.
%
%   Examples:
%
%   Determine the quaternion from rotation angles:
%      yaw = 0.7854; 
%      pitch = 0.1; 
%      roll = 0;
%      q = angle2quat( yaw, pitch, roll )
%
%   Determine the quaternions from multiple rotation angles:
%      yaw = [0.7854 0.5]; 
%      pitch = [0.1 0.3]; 
%      roll = [0 0.1];
%      q = angle2quat( pitch, roll, yaw, 'YXZ' )
%
%   See also DCM2QUAT, QUAT2DCM, QUAT2ANGLE.

%   Copyright 2000-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2010/09/28 03:13:34 $

error(nargchk(3, 4, nargin,'struct'));

if any(~isreal(r1) || ~isnumeric(r1))
    error(message('aero:angle2quat:isNotReal1'));
end

if any(~isreal(r2) || ~isnumeric(r2))
    error(message('aero:angle2quat:isNotReal2'));
end

if any(~isreal(r3) || ~isnumeric(r3))
    error(message('aero:angle2quat:isNotReal3'));
end

if (length(r1) ~= length(r2)) || (length(r1) ~= length(r3))
    error(message('aero:angle2quat:wrongDimension'));
end

if nargin == 3
    type = 'zyx';
else
    if ischar( varargin{1} )
        type = varargin{1};
    else
        error(message('aero:angle2quat:notChar'));
    end
end

angles = [r1(:) r2(:) r3(:)];

cang = cos( angles/2 );
sang = sin( angles/2 );

switch lower( type )
    case 'zyx'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*cang(:,2).*sang(:,3), ...
            sang(:,1).*cang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3)];
    case 'zyz'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            sang(:,1).*cang(:,2).*cang(:,3) + cang(:,1).*cang(:,2).*sang(:,3)];
    case 'zxy'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3)];
    case 'zxz'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            sang(:,1).*sang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3)];
    case 'yxz'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*cang(:,2).*sang(:,3), ...
            sang(:,1).*cang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3)];
    case 'yxy'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            sang(:,1).*cang(:,2).*cang(:,3) + cang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3)];
    case 'yzx'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3)];
    case 'yzy'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            sang(:,1).*sang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3)];
    case 'xyz'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*sang(:,2).*cang(:,3)];
    case 'xyx'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            sang(:,1).*sang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3)];
    case 'xzy'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            sang(:,1).*cang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*cang(:,2).*sang(:,3)];
    case 'xzx'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3)];
    otherwise
        error(message('aero:angle2quat:unknownRotation', type));
end



function qout = quatmultiply( q, varargin )
%  QUATMULTIPLY Calculate the product of two quaternions.
%   N = QUATMULTIPLY( Q, R ) calculates the quaternion product, N, for two
%   given quaternions, Q and R.  Inputs Q and R can be either M-by-4 matrices 
%   containing M quaternions, or a single 1-by-4 quaternion.  N returns an 
%   M-by-4 matrix of quaternion products.  Each element of Q and R must be a
%   real number.  Additionally, Q and R have their scalar number as the first 
%   column.
%
%   Examples:
%
%   Determine the product of two 1-by-4 quaternions:
%      q = [1 0 1 0];
%      r = [1 0.5 0.5 0.75];
%      mult = quatmultiply(q, r)
%
%   Determine the product of a 1-by-4 quaternion with itself:
%      q = [1 0 1 0];
%      mult = quatmultiply(q)
%
%   Determine the product of 1-by-4 and 2-by-4 quaternions:
%      q = [1 0 1 0];
%      r = [1 0.5 0.5 0.75; 2 1 0.1 0.1];
%      mult = quatmultiply(q, r)
%
%   See also QUATCONJ, QUATDIVIDE, QUATINV, QUATMOD, QUATNORM, 
%   QUATNORMALIZE, QUATROTATE.

%   Copyright 2000-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2010/09/28 03:14:27 $

%   Note: Quaternion multiplication is not commutative.

error(nargchk(1, 2, nargin,'struct'));

if any(~isreal(q(:)))
    error(message('aero:quatnorm:isNotReal1'));
end

if (size(q,2) ~= 4)
    error(message('aero:quatnorm:wrongDimension1'));
end

if nargin == 1
    r = q;
else
    r = varargin{1};
    if any(~isreal(r(:)))
        error(message('aero:quatnorm:isNotReal2'));
    end
    if (size(r,2) ~= 4)
        error(message('aero:quatnorm:wrongDimension2'));
    end
    if (size(r,1) ~= size(q,1) && ~( size(r,1) == 1 || size(q,1) == 1))
         error(message('aero:quatnorm:wrongDimension3'));
    end
end

% Calculate vector portion of quaternion product
% vec = s1*v2 + s2*v1 + cross(v1,v2)
vec = [q(:,1).*r(:,2) q(:,1).*r(:,3) q(:,1).*r(:,4)] + ...
         [r(:,1).*q(:,2) r(:,1).*q(:,3) r(:,1).*q(:,4)]+...
         [ q(:,3).*r(:,4)-q(:,4).*r(:,3) ...
           q(:,4).*r(:,2)-q(:,2).*r(:,4) ...
           q(:,2).*r(:,3)-q(:,3).*r(:,2)];

% Calculate scalar portion of quaternion product
% scalar = s1*s2 - dot(v1,v2)
scalar = q(:,1).*r(:,1) - q(:,2).*r(:,2) - ...
             q(:,3).*r(:,3) - q(:,4).*r(:,4);
    
qout = [scalar  vec];



function qout = quatconj( qin ) 
%  QUATCONJ Calculate the conjugate of a quaternion.
%   N = QUATCONJ( Q ) calculates the conjugate, N, for a given quaternion, Q.  
%   Input Q is an M-by-4 matrix containing M quaternions.  N returns an 
%   M-by-4 matrix of conjugates.  Each element of Q must be a real number.  
%   Additionally, Q has its scalar number as the first column.
%
%   Examples:
%
%   Determine the conjugate of q = [1 0 1 0]:
%      conj = quatconj([1 0 1 0])
%
%   See also QUATDIVIDE, QUATINV, QUATMOD, QUATMULTIPLY, QUATNORM, 
%   QUATNORMALIZE, QUATROTATE.

%   Copyright 2000-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2010/09/28 03:14:26 $

if any(~isreal(qin(:)))
    error(message('aero:quatconj:isNotReal'));
end

if (size(qin,2) ~= 4)
    error(message('aero:quatconj:wrongDimension'));
end

qout = [ qin(:,1)  -qin(:,2:4) ];



function varargout=scalar_field_gradient(F,varargin);
% This function returns a numerical approximation to the gradient of the
% N-D size function array 'F' using a variety of explicit finite different
% methods.  The input arguments may be in several forms, however the ouput
% arguments must always be formatted the same.  An output array must be
% specified for each dimension in 'F'.  The output arrays will each be the
% same size as 'F'.  The first input argument form is
%
%  [dFdX1,dFdX2,...,dFdXN] = scalar_field_gradient(F,METHOD)
%
% where 'F' is the N-D array and 'dFdX1', 'dFdX2', et cetera correspond to
% the numerical gradient approximations in the 1st, 2nd, and so on
% dimensions.  In this form the spacing (id est the differences in the
% coordinates) is uniformly set equal to one.  The argument 'METHOD' is
% optional.  If this argument is not specified, the gradient calculation
% method is set to the default value of '2nd_order_central'.
%
% The second input argument form is
%
%  [dFdX1,dFdX2,...,dFdXN] = scalar_field_gradient(F,dX1,dX2,...,dXN,METHOD)
%
% where 'dX1', 'dX2', et cetera are scalar values giving the coordinate
% spacing in the 1st, 2nd, et cetera dimensions.  The argument 'METHOD' is
% optional.  If this argument is not specified, the gradient calculation
% method is set to the default value of '2nd_order_central'.
%
% The third input argument form is
%
%  [dFdX1,dFdX2,...,dFdXN] = scalar_field_gradient(F,X1,X2,...,XN,METHOD)
%
% where 'X1', 'X2', et cetera are vectors giving the coordinates of the
% function values in 'F' in each dimension.  The lengths of the vectors
% 'X1', 'X2', et cetera must equal the size of 'F' in the corresponding
% dimension.  Id est
%
%  length(X1) = size(F,1)
%  length(X2) = size(F,2)
%   ...
%  length(XN) = size(F,N)
%
% The argument 'METHOD' is optional.  If this argument is not specified,
% the gradient calculation method is set to the default value of
% '2nd_order_central'.
%
% The fourth input argument form is
%
%  [dFdX1,dFdX2,...,dFdXN] = scalar_field_gradient(F,XA1,XA2,...,XAN,METHOD)
%
% where 'X1', 'X2', et cetera are N-D arrays of the same size as 'F' which
% give the coordinates of the function values in 'F' in each dimension.
% The argument 'METHOD' is optional.  If this argument is not specified,
% the gradient calculation method is set to the default value of
% '2nd_order_central'.
%
% The gradient calculation arguent 'METHOD' may be set equal to one of the
% following strings
%
%  '2nd_order_central'                  0.710*sigma/dx
%  '4th_order_central'                  0.950*sigma/dx
%  '2nd_order_richardson_noise_min'     0.088*sigma/dx
%  '4th_order_richardson_noise_min'     0.334*sigma/dx
%  '2nd_order_least_squares'            0.316*sigma/dx
%  '6th_order_implicit'                 1.000*sigma/dx
%  '10th_order_implicit'                0.840*sigma/dx
%
% where the second column gives the expected noise error in terms of the
% measurement uncertainty 'sigma' and the grid spacing 'dx'.  These methods
% were taken from section 9.3 of 'Particle Image Velocimetry' by Ronald J.
% Adrian and Jerry Westerweel and are discussed more fully in the paper
% 'Some considerations on the accuracy and frequency response of some
% derivative filters applied to particle image velocimetry fields' by J. M.
% Foucaut and M. Stanislas 2002.
%
% The derivative estimates calculated near the edge of the array 'F' where
% there are not a sufficient number of function values to use the specified
% central difference methods use standard forward/backward finite
% difference methods.  If the array 'F' is too small to use a specified
% explicit finite difference method in any dimension, then the appropriate
% order forward/backward difference methods will be used in this dimension.
% If the array 'F' is too small to use a specified implicit finite
% difference method in any dimension, then the 2nd order central difference
% will be used for central values and the appropriate forward/backward
% difference will be used for components near the edge of the data field.
%
% Author: Rod La Foy
% First Written On: 22 April 2014
% Last Modified On: 25 April 2014

% This is the size of the input array 'F'
array_size=size(F);
% This is the number of dimensions of the input array 'F'
array_dimension_number=ndims(F);

% This checks whether the input array is 2D
if array_dimension_number==2;
    % This checks whether the input array is actually 1D and then removes
    % the singular dimension
    if array_size(2)==1;
        % This removes the second dimension size
        array_size(2)=[];
        % This sets the number of dimensions in the intput array to 1
        array_dimension_number=1;
    end;
end;

% This is the number of elements input into the variable length input
% argument list 'varargin'
argument_number=length(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This checks that the input/output arguments are formatted correctly.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If the number of input arguments in 'varargin' is less then the number of
% dimensions of the input function array 'F' then an error is returned
if ((argument_number<array_dimension_number)&&(argument_number~=0))&&not((argument_number==1)&&(ischar(varargin{1})));
    % This returns an error stating that an incorrect number of inputs was
    % entered
    error(['An incorrect number of arguments was input.  The input array ''F'' has ',num2str(array_dimension_number),' dimensions, however only ',num2str(argument_number+1),' arguments were input.']);
end;

% If the number of output arguments 'nargout' does not equal the number of
% dimension in the input array 'F' then an error is returned
if nargout~=array_dimension_number;
    % This returns an error stating that an incorrect number of output
    % arguments were used
    error(['An incorrect number of ouput arguments was used.  The input array ''F'' has ',num2str(array_dimension_number),' dimensions which should equal the number of output arguments.']);
end;

% If the number of input arguments in 'varargin' equals the number of
% dimensions of the input array 'F' then the calculation method is assummed
% to be equal to the default method
if argument_number==array_dimension_number;
    % This sets the calculation method equal to the default method
    gradient_method='default';
end;

% If there is only one optional argument and it is a string, then this sets
% the gradient calculation method equal to the string
if (argument_number==1)&&(ischar(varargin{1}));
    % This sets the gradient calculation method equal to the last
    % member of the cell 'varargin'
    gradient_method=varargin{1};
    % This creates a Boolean value to check whether the input string
    % defining the gradient calculation method matches one of the
    % defined methods
    gradient_method_check=not(strcmp(gradient_method,'2nd_order_central'));
    gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'4th_order_central'));
    gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'2nd_order_richardson_noise_min'));
    gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'4th_order_richardson_noise_min'));
    gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'2nd_order_least_squares'));
    gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'6th_order_implicit'));
    gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'10th_order_implicit'));
    % This checks that the gradient method equals one of the defined
    % calculation methods
    if gradient_method_check;
        % This displays an error stating that the gradient calculation
        % method string does not equal one of the defined methods
        error('The first optional input argument does not equal one of the defined gradient calculation methods.');
    end;
end;

% If the number of input arguments in 'varargin' equals the one more than
% the number of dimensions of the input array 'F' then the calculation
% method is set equal to the last input argument (if it is a string)
if argument_number==(array_dimension_number+1);
    % This checks whether the last member of the cell array 'varargin' is a
    % string and if so, set the input method equal to the string and if
    % not, returns an error
    if ischar(varargin{end});
        % This sets the gradient calculation method equal to the last
        % member of the cell 'varargin'
        gradient_method=varargin{end};
        % This creates a Boolean value to check whether the input string
        % defining the gradient calculation method matches one of the
        % defined methods
        gradient_method_check=not(strcmp(gradient_method,'2nd_order_central'));
        gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'4th_order_central'));
        gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'2nd_order_richardson_noise_min'));
        gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'4th_order_richardson_noise_min'));
        gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'2nd_order_least_squares'));
        gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'6th_order_implicit'));
        gradient_method_check=gradient_method_check&&not(strcmp(gradient_method,'10th_order_implicit'));
        % This checks that the gradient method equals one of the defined
        % calculation methods
        if gradient_method_check;
            % This displays an error stating that the gradient calculation
            % method string does not equal one of the defined methods
            error('The last input argument does not equal one of the defined gradient calculation methods.');
        end;
    else;
        % This displays an error stating that the last element of the
        % 'varargin' cell array is not a string and thus cannot be giving
        % the calculation method
        error('The last input argument was not a string.  Either there are an incorrect number of inputs or the input arguments are not formatted correctly.');
    end;
end;

% If the number of input arguments in 'varargin' is more than one more than
% the number of dimensions of the input array 'F' then an error is returned
if argument_number>(array_dimension_number+1);
    % This returns an error stating that an incorrect number of inputs was
    % entered
    error(['An incorrect number of arguments was input.  The input array ''F'' has ',num2str(array_dimension_number),' dimensions, thus at most ',num2str(array_dimension_number+2),' arguments can be input.']);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This extracts the grid spacing from the input arguments.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This initializes a variable to measure the size of the input independent
% variables
independent_variable_size=zeros(array_dimension_number,array_dimension_number);

% If the only input argument was the array 'F' then the spacing is set
% equal to 1, otherwise, this calculates the array spacing based upon the
% input arguments
if (argument_number==0)||((argument_number==1)&&(ischar(varargin{1})));
    
    % This sets the calculation method equal to the default method
    gradient_method='default';
    
    % This sets the spacing values of the independent variables to all
    % equal one
    dX=ones(array_dimension_number,1);
    
else;
    
    % This iterates through the dimensions of 'F' measuring the size of the
    % indepdent variables to determine what for of input the data is in
    for dimension_index=1:array_dimension_number;
        % This is the size of the current independent variable
        current_independent_variable_size=size(varargin{dimension_index});
        % This is the size of the current independent variable
        independent_variable_size(dimension_index,1:length(current_independent_variable_size))=current_independent_variable_size;
    end;
    
    % This initializes the vector of spacing values of the independent
    % variables
    dX=zeros(array_dimension_number,1);
    
    % This creates a vector equal to [1,1,0,0,...] which has the same number of
    % elements as the number of dimensions in the array 'F' to test whether the
    % input arguments in 'varargin' are all scalars
    scalar_test_vector=zeros(1,array_dimension_number);
    % This sets the first two elements equal to 1
    scalar_test_vector(1:2)=1;
    
    % This tests whether all the dimensions of the independent variables are
    % equal to [1,1] which would imply that the input variables are scalars
    % equal to the spacing in each dimension
    if all(all(bsxfun(@eq,independent_variable_size,scalar_test_vector)));
        
        % This iterates through the input arguments of 'varargin' extracting
        % the scalar spacing values
        for dimension_index=1:array_dimension_number;
            % This is the current independent variable spacing
            dX(dimension_index)=varargin{dimension_index};
        end;
        
        % This tests whether all the independent variables are vectors which would
        % imply that the spacing can be taken as the difference between the
        % elements first and second elements of the vector
    elseif all(sum(bsxfun(@gt,independent_variable_size,1),2)==ones(array_dimension_number,1));
        
        % This iterates through the input arguments of 'varargin' calculating
        % the spacing values
        for dimension_index=1:array_dimension_number;
            
            % This is the current independent variable
            X_Current=varargin{dimension_index};
            
            % This tests whether the length of the vector equals the number
            % size of the input array in the current dimension and if not
            % displays an error
            if length(X_Current)~=array_size(dimension_index);
                % This displays an error stating that the length of the current
                % independent vector does not equal the size of the array 'F'
                % in the current dimension
                error('A coordinate vector does not have the same length as the size of the array ''F'' in at least one dimension.');
            end;
            
            % This is the spacing between the elements of the vector
            dX_Current=diff(X_Current);
            % This calculates the number of decimal points of precision
            % used to specify the spacing
            decimal_point_number=floor(log10(dX_Current(1)/eps(dX_Current(1))));
            % This rounds the spacing to the the machine precision
            dX_Current=precision_round(dX_Current,decimal_point_number-3);
            
            % This tests whether there is more than one value of the spacings
            % and returns an error if so
            if length(unique(dX_Current))>1;
                % This displays an error stating that the spacing is not
                % constant in the current dimension
                error('A coordinate vector does not have constant spacing between it''s elements.');
            end;
            
            % This is the current independent variable spacing
            dX(dimension_index)=dX_Current(1);
            
        end;
        
        % This clears the temporary arrays from memory
        clear('X_Current','dX_Current');
        
        % This tests whether all the independent variables are arrays of the same
        % size as the input array 'F' and then sets the spacing as the difference
        % between the first and second elements in each dimension
    elseif all(all(bsxfun(@eq,independent_variable_size,array_size)));
        
        % This iterates through the input arguments of 'varargin' calculating
        % the spacing values
        for dimension_index=1:array_dimension_number;
            
            % This is the current independent variable
            X_Current=varargin{dimension_index};
            
            % This is the spacing between the elements of the array in the
            % current dimension
            dX_Current=diff(X_Current,1,dimension_index);
            % This calculates the number of decimal points of precision
            % used to specify the spacing
            decimal_point_number=floor(log10(dX_Current(1)/eps(dX_Current(1))));
            % This rounds the spacing to the the machine precision
            dX_Current=precision_round(dX_Current,decimal_point_number-3);
            
            % This tests whether there is more than one value of the spacings
            % and returns an error if so
            if length(unique(dX_Current))>1;
                
                % This displays an error stating that the spacing is not
                % constant in the current dimension
                error('A coordinate array does not have constant spacing between it''s elements.');
            end;
            
            % This is the current independent variable spacing
            dX(dimension_index)=dX_Current(1);
            
        end;
        
        % This clears the temporary arrays from memory
        clear('X_Current','dX_Current');
        
    else;
        
        % This displays an error stating input independent variables are not
        % formatted correctly
        error('The input coordinates are not formatted correctly.');
        
    end;
    
    % This tests whether any of the spacing variables are zero and if so
    % returns an error
    if any(dX==0);
        % This displays an error stating that the some of the spacing variables
        % are zero
        error('At least one spacing variable ''dX'' is equal to zero.  This may occur if ''meshgrid'' was used to generate the coordinates rather than ''ndgrid''.  This may be fixed by permuting the 1st and 2nd dimensions.');
    end;
    
end;

% This initializes the output cell array
varargout=cell(array_dimension_number,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This specifies the finite difference kernals to use.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This calculates the first derivative using the second order central
% difference
if strcmp(gradient_method,'2nd_order_central')||strcmp(gradient_method,'default');
    
    % This is the 2nd order central difference coefficient kernal
    center_coefficient_kernal=[-1/2,0,+1/2];
    % This is the 2nd order central difference coefficient relative
    % coordinates
    center_coefficient_coordinates=[-1,0,+1];
    
    % This is the order of the finite difference estimate
    finite_difference_order=2;
    
    % This sets the computational method as explicit
    computational_method='explicit';
    
    % This calculates the first derivative using the fourth order central
    % difference
elseif strcmp(gradient_method,'4th_order_central');
    
    % This is the 4th order central difference coefficient kernal
    center_coefficient_kernal=[+1/12,-2/3,0,+2/3,-1/12];
    % This is the 4th order central difference coefficient relative
    % coordinates
    center_coefficient_coordinates=[-2,-1,0,+1,+2];
    
    % This is the order of the finite difference estimate
    finite_difference_order=4;
    
    % This sets the computational method as explicit
    computational_method='explicit';
    
    % This calculates the first derivative using the second order Richardson
    % extrapolation
elseif strcmp(gradient_method,'2nd_order_richardson_noise_min');
    
    % This is the 2nd order Richardson difference coefficient kernal
    center_coefficient_kernal=[-4/65,0,0,0,0,0,0,+1/130,0,-1/130,0,0,0,0,0,0,+4/65];
    % This is the 2nd order Richardson difference coefficient relative
    % coordinates
    center_coefficient_coordinates=[-8,-7,-6,-5,-4,-3,-2,-1,0,+1,+2,+3,+4,+5,+6,+7,+8];
    
    % This is the order of the finite difference estimate
    finite_difference_order=2;
    
    % This sets the computational method as explicit
    computational_method='explicit';
    
    % This calculates the first derivative using the fourth order Richardson
    % extrapolation
elseif strcmp(gradient_method,'4th_order_richardson_noise_min');
    
    % This is the 4th order Richardson difference coefficient kernal
    center_coefficient_kernal=[+23/6608,0,0,0,0,0,-37/177,-136/1239,0,+136/1239,+37/177,0,0,0,0,0,-23/6608];
    % This is the 4th order Richardson difference coefficient relative
    % coordinates
    center_coefficient_coordinates=[-8,-7,-6,-5,-4,-3,-2,-1,0,+1,+2,+3,+4,+5,+6,+7,+8];
    
    % This is the order of the finite difference estimate
    finite_difference_order=4;
    
    % This sets the computational method as explicit
    computational_method='explicit';
    
    % This calculates the first derivative using the second order least
    % squares central difference
elseif strcmp(gradient_method,'2nd_order_least_squares');
    
    % This is the 2nd order central difference coefficient kernal
    center_coefficient_kernal=[-1/5,-1/10,0,+1/10,+1/5];
    % This is the 2nd order central difference coefficient relative
    % coordinates
    center_coefficient_coordinates=[-2,-1,0,+1,+2];
    
    % This is the order of the finite difference estimate
    finite_difference_order=2;
    
    % This sets the computational method as explicit
    computational_method='explicit';
    
    % This calculates the first derivative using the sixth order implicit
    % compact difference method
elseif strcmp(gradient_method,'6th_order_implicit');
    
    % This is the derivative coefficient kernal for the 6th order
    % implicit method
    derivative_center_coefficient_kernal=[+1/3,+1,+1/3];
    % This is the function coefficient kernal for the 6th order implicit
    % method
    function_center_coefficient_kernal=[-1/36,-7/9,0,+7/9,+1/36];
    
    % This is the 2nd order central difference coefficient kernal
    center_coefficient_kernal=[-1/2,0,+1/2];
    % This is the 2nd order central difference coefficient relative
    % coordinates
    center_coefficient_coordinates=[-1,0,+1];
    
    % This is the order of the finite difference estimate
    finite_difference_order=2;
    
    % This sets the computational method as implicit
    computational_method='implicit';
    
    % This calculates the first derivative using the tenth order implicit
    % compact difference method
elseif strcmp(gradient_method,'10th_order_implicit');
    
    % This is the derivative coefficient kernal for the 10th order
    % implicit method
    derivative_center_coefficient_kernal=[+1/20,+1/2,+1,+1/2,+1/20];
    % This is the function coefficient kernal for the 10th order implicit
    % method
    function_center_coefficient_kernal=[-1/600,-101/600,-17/24,0,17/24,+101/600,+1/600];
    
    % This is the 2nd order central difference coefficient kernal
    center_coefficient_kernal=[-1/2,0,+1/2];
    % This is the 2nd order central difference coefficient relative
    % coordinates
    center_coefficient_coordinates=[-1,0,+1];
    
    % This is the order of the finite difference estimate
    finite_difference_order=2;
    
    % This sets the computational method as implicit
    computational_method='implicit';
    
end;

% This creates a table of coefficients for the forward/backward finite
% difference estimates
coefficient_table_values{1}=[  -1,         1       ];
coefficient_table_values{2}=[  -3/2,       2,      -1/2        ];
coefficient_table_values{3}=[  -11/6,      3,      -3/2,       1/3       ];
coefficient_table_values{4}=[  -25/12,     4,      -3,         4/3,        -1/4       ];

% This creates a table of the coordinates of the coefficients of the
% forward/backward finite difference estimates
coefficient_table_coordinates{1}=0:1;
coefficient_table_coordinates{2}=0:2;
coefficient_table_coordinates{3}=0:3;
coefficient_table_coordinates{4}=0:4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This calculates the finite difference estimate in each dimension.       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This iterates through the different dimension calculating the
% derivatives using the defined finite difference scheme
for dimension_index=1:array_dimension_number;
    
    % This calculates the finite difference estimate using an explicit
    % method
    if strcmp(computational_method,'explicit');
        
        % This is the domain in the current dimension over which the
        % central difference method may be used
        central_difference_domain_min=length(center_coefficient_kernal);
        central_difference_domain_max=array_size(dimension_index)-(length(center_coefficient_kernal)-1);
        
        % If the domain does not fit within the array volume, then the
        % values of the domain minimum and maximum are set such that volume
        % is split into calculating the derivative of half of the array
        % with forward difference and half with backward difference
        if (central_difference_domain_min>array_size(dimension_index))||(central_difference_domain_max<1);
            % This sets the domain such that the volume is split in the
            % middle
            central_difference_domain_min=round(array_size(dimension_index)/2)+1;
            central_difference_domain_max=round(array_size(dimension_index)/2);
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This calculates the central difference data coefficients.           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % This is the domain over which the finite difference may be calculated
        calculation_index_min=ones(1,array_dimension_number);
        calculation_index_max=array_size;
        
        % This sets the domain to minimum and maximum domain indices in
        % the current dimension
        calculation_index_min(dimension_index)=central_difference_domain_min;
        calculation_index_max(dimension_index)=central_difference_domain_max;
        
        % This is the current coefficient kernal to use for calculating the
        % finite difference estimate
        scaled_coefficient_kernal=center_coefficient_kernal/dX(dimension_index);
        
        % This calculates the finite difference coefficient sparse array
        finite_difference_coefficients=calculate_coefficient_array(array_size,calculation_index_min,calculation_index_max,dimension_index,scaled_coefficient_kernal,center_coefficient_coordinates);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This calculates the forward difference data coefficients.           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % This is the domain over which the forward difference may be
        % calculated
        calculation_index_min=ones(1,array_dimension_number);
        calculation_index_max=array_size;
        
        % This iterates through the range of indices that need to be
        % calculated using the forward difference method
        for calculation_index_temp=1:(central_difference_domain_min-1);
            
            % This sets the domain to current index in the current
            % dimension
            calculation_index_min(dimension_index)=calculation_index_temp;
            calculation_index_max(dimension_index)=calculation_index_temp;
            
            % This is the number of elements available to calculate the finite
            % difference estimate
            available_element_number=array_size(dimension_index)-calculation_index_max(dimension_index)+1;
            
            % This is the order of the coefficient kernal to use in the
            % forward difference estimate
            current_finite_difference_order=min([finite_difference_order,available_element_number-1]);
            
            % This is the current coefficient kernal to use for calculating the
            % forward finite difference estimate
            scaled_coefficient_kernal=coefficient_table_values{current_finite_difference_order}/dX(dimension_index);
            
            % These are the forward difference coefficient coordinates
            forward_coefficient_coordinates=coefficient_table_coordinates{current_finite_difference_order};
            
            % This calculates the forward difference coefficient sparse array
            finite_difference_coefficients=finite_difference_coefficients+calculate_coefficient_array(array_size,calculation_index_min,calculation_index_max,dimension_index,scaled_coefficient_kernal,forward_coefficient_coordinates);
            
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This calculates the backward difference data coefficients.          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % This is the domain over which the backward difference may be
        % calculated
        calculation_index_min=ones(1,array_dimension_number);
        calculation_index_max=array_size;
        
        % This iterates through the range of indices that need to be
        % calculated using the backward difference method
        for calculation_index_temp=(central_difference_domain_max+1):array_size(dimension_index);
            
            % This sets the domain to current index in the current
            % dimension
            calculation_index_min(dimension_index)=calculation_index_temp;
            calculation_index_max(dimension_index)=calculation_index_temp;
            
            % This is the number of elements available to calculate the finite
            % difference estimate
            available_element_number=calculation_index_min(dimension_index);
            
            % This is the order of the coefficient kernal to use in the
            % backward difference estimate
            current_finite_difference_order=min([finite_difference_order,available_element_number-1]);
            
            % This is the current coefficient kernal to use for calculating the
            % backward finite difference estimate
            scaled_coefficient_kernal=-coefficient_table_values{current_finite_difference_order}/dX(dimension_index);
            
            % These are the forward difference coefficient coordinates
            backward_coefficient_coordinates=-coefficient_table_coordinates{current_finite_difference_order};
            
            % This calculates the backward difference coefficient sparse
            % array
            finite_difference_coefficients=finite_difference_coefficients+calculate_coefficient_array(array_size,calculation_index_min,calculation_index_max,dimension_index,scaled_coefficient_kernal,backward_coefficient_coordinates);
            
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This calculates the full finite difference.                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % This calculates the finite difference estimate in the current
        % dimension
        varargout{dimension_index}=reshape(finite_difference_coefficients*F(:),array_size);
        
        % This clears the current coefficient sparse array from memory
        clear('finite_difference_coefficients');
        
    % This calculates the finite difference estimate using an implicit
    % method    
    elseif strcmp(computational_method,'implicit');
        
        % This is the "radius" of the function coefficient kernal
        function_coefficient_radius=((length(function_center_coefficient_kernal)-1)/2)+1;
        
        % This is the domain in the current dimension over which the
        % central difference method may be used
        central_difference_domain_min=max([length(derivative_center_coefficient_kernal),length(function_center_coefficient_kernal),length(coefficient_table_values{finite_difference_order})]);
        central_difference_domain_max=array_size(dimension_index)-(max([length(derivative_center_coefficient_kernal),length(function_center_coefficient_kernal),length(coefficient_table_values{finite_difference_order})])-1);
        
        % If the domain does not fit within the array volume, then the
        % values of the domain minimum and maximum are set such that volume
        % is split into calculating the derivative of half of the array
        % with forward difference and half with backward difference
        if (central_difference_domain_min>array_size(dimension_index))||(central_difference_domain_max<1);
            % This sets the domain such that the volume is split in the
            % middle
            central_difference_domain_min=round(array_size(dimension_index)/2)+1;
            central_difference_domain_max=round(array_size(dimension_index)/2);
        end;
        
        % This checks that there are enough components of the array in the
        % current dimension to use the implicit calculation scheme
        if central_difference_domain_min<central_difference_domain_max;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This calculates the implicit finite difference.             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % This creates the diagonal sparse array for the derivative coefficients
            derivative_coefficient_diagonal_array=create_diagonal_matrix(array_size(dimension_index),derivative_center_coefficient_kernal);
            
            % This creates the diagonal sparse array for the function coefficients
            function_coefficient_diagonal_array=create_diagonal_matrix(array_size(dimension_index),function_center_coefficient_kernal);
            
            % This iterates through the sparse array adding in coefficients for solving
            % for the derivatives of the edge terms using an explicit forward
            % difference method
            for ii_coefficient_index=1:(function_coefficient_radius-1);
                
                % This initializes a vector to substitute into the derivative
                % coefficient array
                derivative_diagonal_array_vector=zeros(1,array_size(dimension_index));
                % This sets the current ii index of the derivative vector to 1 (ie
                % create a vector for explicity calculating the edge component
                % derivatives)
                derivative_diagonal_array_vector(ii_coefficient_index)=1;
                % This sets the current ii index row vector of the derivative sparse
                % array equal to the defined vector
                derivative_coefficient_diagonal_array(ii_coefficient_index,:)=derivative_diagonal_array_vector;
                
                % This initializes a vector to substitute into the function
                % coefficient array
                function_diagonal_array_vector=zeros(1,array_size(dimension_index));
                % This sets the function vector equal to the coefficients for
                % calculating a forward finite difference
                function_diagonal_array_vector(ii_coefficient_index+(0:(length(coefficient_table_values{finite_difference_order})-1)))=coefficient_table_values{finite_difference_order};
                % This sets the current ii index row vector of the function sparse
                % array equal to the defined vector
                function_coefficient_diagonal_array(ii_coefficient_index,:)=function_diagonal_array_vector;
                
            end;
            
            % This iterates through the sparse array adding in coefficients for solving
            % for the derivatives of the edge terms using an explicit backward
            % difference method
            for ii_coefficient_index=(array_size(dimension_index)-(function_coefficient_radius-2)):array_size(dimension_index);
                
                % This initializes a vector to substitute into the derivative
                % coefficient array
                derivative_diagonal_array_vector=zeros(1,array_size(dimension_index));
                % This sets the current ii index of the derivative vector to 1 (ie
                % create a vector for explicity calculating the edge component
                % derivatives)
                derivative_diagonal_array_vector(ii_coefficient_index)=1;
                % This sets the current ii index row vector of the derivative sparse
                % array equal to the defined vector
                derivative_coefficient_diagonal_array(ii_coefficient_index,:)=derivative_diagonal_array_vector;
                
                % This initializes a vector to substitute into the function
                % coefficient array
                function_diagonal_array_vector=zeros(1,array_size(dimension_index));
                % This sets the function vector equal to the coefficients for
                % calculating a backward finite difference
                function_diagonal_array_vector(ii_coefficient_index-(0:(length(coefficient_table_values{finite_difference_order})-1)))=-coefficient_table_values{finite_difference_order};
                % This sets the current ii index row vector of the function sparse
                % array equal to the defined vector
                function_coefficient_diagonal_array(ii_coefficient_index,:)=function_diagonal_array_vector;
                
            end;
            
            % This calculates the full implicit derivative coefficient array
            coefficient_array=full(inv(derivative_coefficient_diagonal_array)*function_coefficient_diagonal_array);
            
            % This is a permutation vector used to permute the function array 'F' so
            % that the current dimension is the first dimension
            permutation_vector=1:array_dimension_number;
            % This sets the first value of the permutation vector to the current
            % dimension
            permutation_vector(1)=dimension_index;
            % This sets the dimension of the permutation vector to 1
            permutation_vector(dimension_index)=1;
            
            % This permutes function array 'F' so that the current processing dimension
            % is in the first dimension
            F_Permute=permute(F,permutation_vector);
            
            % This initializes the current array of the output cell array
            varargout{dimension_index}=zeros(array_size(permutation_vector));
            
            % This iterates through all the other dimensions of the function array 'F'
            % calculating the derivate of the current vector
            for vector_index=1:prod([array_size(1:(dimension_index-1)),array_size((dimension_index+1):end)]);
                
                % This is the domain of the current vector to exract from the function
                % array 'F'
                ii_index_min=1+(vector_index-1)*array_size(dimension_index);
                ii_index_max=vector_index*array_size(dimension_index);
                
                % This extracts the current vector from the function array 'F'
                F_Vector=F_Permute(ii_index_min:ii_index_max)';
                
                % This calculates the derivative of the current column of the function
                % array 'F' using the implicit coefficient array
                varargout{dimension_index}(ii_index_min:ii_index_max)=(coefficient_array*F_Vector)*(1/dX(dimension_index));
                
            end;
            
            % This inverts the permutation of the output derivative array
            varargout{dimension_index}=ipermute(varargout{dimension_index},permutation_vector);
            
        else;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This calculates the central difference data coefficients.   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % This is the domain over which the finite difference may be calculated
            calculation_index_min=ones(1,array_dimension_number);
            calculation_index_max=array_size;
            
            % This sets the domain to minimum and maximum domain indices in
            % the current dimension
            calculation_index_min(dimension_index)=central_difference_domain_min;
            calculation_index_max(dimension_index)=central_difference_domain_max;
            
            % This is the current coefficient kernal to use for calculating the
            % finite difference estimate
            scaled_coefficient_kernal=center_coefficient_kernal/dX(dimension_index);
            
            % This calculates the finite difference coefficient sparse array
            finite_difference_coefficients=calculate_coefficient_array(array_size,calculation_index_min,calculation_index_max,dimension_index,scaled_coefficient_kernal,center_coefficient_coordinates);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This calculates the forward difference data coefficients.   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % This is the domain over which the forward difference may be
            % calculated
            calculation_index_min=ones(1,array_dimension_number);
            calculation_index_max=array_size;
            
            % This iterates through the range of indices that need to be
            % calculated using the forward difference method
            for calculation_index_temp=1:(central_difference_domain_min-1);
                
                % This sets the domain to current index in the current
                % dimension
                calculation_index_min(dimension_index)=calculation_index_temp;
                calculation_index_max(dimension_index)=calculation_index_temp;
                
                % This is the number of elements available to calculate the finite
                % difference estimate
                available_element_number=array_size(dimension_index)-calculation_index_max(dimension_index)+1;
                
                % This is the order of the coefficient kernal to use in the
                % forward difference estimate
                current_finite_difference_order=min([finite_difference_order,available_element_number-1]);
                
                % This is the current coefficient kernal to use for calculating the
                % forward finite difference estimate
                scaled_coefficient_kernal=coefficient_table_values{current_finite_difference_order}/dX(dimension_index);
                
                % These are the forward difference coefficient coordinates
                forward_coefficient_coordinates=coefficient_table_coordinates{current_finite_difference_order};
                
                % This calculates the forward difference coefficient sparse array
                finite_difference_coefficients=finite_difference_coefficients+calculate_coefficient_array(array_size,calculation_index_min,calculation_index_max,dimension_index,scaled_coefficient_kernal,forward_coefficient_coordinates);
                
            end;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This calculates the backward difference data coefficients.  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % This is the domain over which the backward difference may be
            % calculated
            calculation_index_min=ones(1,array_dimension_number);
            calculation_index_max=array_size;
            
            % This iterates through the range of indices that need to be
            % calculated using the backward difference method
            for calculation_index_temp=(central_difference_domain_max+1):array_size(dimension_index);
                
                % This sets the domain to current index in the current
                % dimension
                calculation_index_min(dimension_index)=calculation_index_temp;
                calculation_index_max(dimension_index)=calculation_index_temp;
                
                % This is the number of elements available to calculate the finite
                % difference estimate
                available_element_number=calculation_index_min(dimension_index);
                
                % This is the order of the coefficient kernal to use in the
                % backward difference estimate
                current_finite_difference_order=min([finite_difference_order,available_element_number-1]);
                
                % This is the current coefficient kernal to use for calculating the
                % backward finite difference estimate
                scaled_coefficient_kernal=-coefficient_table_values{current_finite_difference_order}/dX(dimension_index);
                
                % These are the forward difference coefficient coordinates
                backward_coefficient_coordinates=-coefficient_table_coordinates{current_finite_difference_order};
                
                % This calculates the backward difference coefficient sparse
                % array
                finite_difference_coefficients=finite_difference_coefficients+calculate_coefficient_array(array_size,calculation_index_min,calculation_index_max,dimension_index,scaled_coefficient_kernal,backward_coefficient_coordinates);
                
            end;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This calculates the full finite difference.                 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%             % This permutes function array 'F' so that the current processing dimension
%             % is in the first dimension
%             F=permute(F,permutation_vector);
            
            
            % This calculates the finite difference estimate in the current
            % dimension
            varargout{dimension_index}=reshape(finite_difference_coefficients*F(:),array_size);
            
            % This clears the current coefficient sparse array from memory
            clear('finite_difference_coefficients');
            
        end;
        
    end;
    
end;



function finite_difference_coefficients=calculate_coefficient_array(array_size,calculation_index_min,calculation_index_max,dimension_index,coefficient_kernal,coefficient_coordinates);
% This function creates the sparse array 'finite_difference_coefficients'
% which when multiplied by an array of size 'array_size' will produce an
% array whose values approximate the derivative of the array.  Specifically
% for a function array 'F', the approximation to the derivative is given by
%
%  dUdX = reshape( finite_difference_coefficients * F(:), array_size );
%
% This input argument 'array_size' is a 1 x N vector giving the size of
% the array to produce the coefficent array for.  The vectors
% 'calculation_index_min' and 'calculation_index_max' are also 1 x N
% vectors that give the minimum and maximum domain over which the
% coefficients are solved.  The argument 'coefficient_kernal' is a 1 x M
% vector giving the weights of the values of the function to the
% derivative.  The argument 'coefficient_coordinates' gives the
% spatial relationship between the coefficient coordinates within the
% function array.

% This is the number of dimensions in the array
array_dimension_number=length(array_size);
% This is the total number of elements in the array
array_element_number=prod(array_size);

% This is the number of coefficients to be used in the calculation
coefficient_number=length(coefficient_kernal);

% This initializes a cell array to store the vectors of the domain
% over which to calculate the 2nd order central difference
coordinate_vector_list=cell(1,array_dimension_number);

% This creates a cell array of vectors corresponding to range of
% indices to calculate the 2nd order central difference over
for vector_index=1:array_dimension_number;
    % This adds the current vector to the cell array
    coordinate_vector_list{vector_index}=calculation_index_min(vector_index):calculation_index_max(vector_index);
end;

% This creates a cell array of coordinate arrays giving the
% subscripted indices of the array at which to calculate the 2nd
% order central difference
coordinate_array_list=ndgrid_vectorized(coordinate_vector_list);

% This is the number of elements in the coordinate arrays
calculation_element_number=numel(coordinate_array_list{1});

% This initializes a subscripted index array giving the subscripted
% coordinates at which to calculate the 2nd order central
% difference
subscript_index_array=zeros(calculation_element_number,array_dimension_number);

% This fills in the subscripted index array with the coordinate
% arrays
for vector_index=1:array_dimension_number;
    % This adds the current dimension's subscripts to the array
    subscript_index_array(:,vector_index)=coordinate_array_list{vector_index}(:);
end;

% This clears the coordinate array list from memory as it may be
% large
clear('coordinate_array_list');

% These initializes the jj subscripted indices of the locations of the
% sparse finite difference coefficient array that are non-zero
coefficient_jj_indices=zeros(calculation_element_number*coefficient_number,1);
% These initializes the finite difference coefficient values of the sparse
% array
coefficient_values=zeros(calculation_element_number*coefficient_number,1);

% This iterates through the different coefficients calculating the linear
% indices
for coefficient_index=1:coefficient_number;
    
    % This copies the subscripted index array into a temporary array for
    % determining the current coefficient
    subscript_index_array_temp=subscript_index_array;
    
    % This offsets the subscripted array by the defined coordinate distance
    % in the current dimension
    subscript_index_array_temp(:,dimension_index)=subscript_index_array_temp(:,dimension_index)+coefficient_coordinates(coefficient_index);
    
    % This calculates the linear indices of the current subscripted index
    % coordinates
    linear_indices_temp=sub2ind_vectorized(array_size,subscript_index_array_temp);
    % This clears the current subscripted index array from memory
    clear('subscript_index_array_temp');
    
    % If the current coefficient corresponds to the 0 coefficient
    % coordinate, then the 'coefficient_ii_indices' are calculated
    if coefficient_coordinates(coefficient_index)==0;
        % This calculates the ii coordinates of the coefficients of the
        % sparse array
        coefficient_ii_indices=repmat(linear_indices_temp,[coefficient_number,1]);
    end;
    
    % This adds the current jj coordinates of the coefficients to the full
    % jj coordinate vector
    coefficient_jj_indices((1:calculation_element_number)+(coefficient_index-1)*calculation_element_number)=linear_indices_temp;
    % This clears the current linear indices from memory
    clear('linear_indices_temp');
    
    % This adds the current coefficient values to the full coefficent value
    % vector
    coefficient_values((1:calculation_element_number)+(coefficient_index-1)*calculation_element_number)=coefficient_kernal(coefficient_index);
    
end;

% This creates the finite difference coefficient array
finite_difference_coefficients=sparse(coefficient_ii_indices,coefficient_jj_indices,coefficient_values,array_element_number,array_element_number);



function gridded_array=ndgrid_vectorized(vector_list);
% This function performs a similar operation to the builtin matlab function
% 'ndgrid' however this function can be applied on an arbitrary number of
% input vectors.  The input vectors are stored within the cell array
% 'vector_list' such that the individual vectors in each dimension are
%
%  x1 = vector_list{1};
%  x2 = vector_list{2};
%  x3 = vector_list{3};
%   ...
%  xN = vector_list{N};
%
% and the output 'gridded_array' is a cell array such that
%
%  gridded_array{1} = X1;
%  gridded_array{2} = X2;
%  gridded_array{3} = X3;
%   ...
%  gridded_array{N} = XN;
%
% where the variables X1, X2, X3, . . . XN would be produced by
%
%  [X1,X2,X3,...,XN] = ndgrid(x1,x2,x3,...,xN)
%
% by the standard 'ndgrid' function.

% This is the number of dimensions to create
array_dimension_number=length(vector_list);

% This initializes the output cell array
gridded_array=cell(array_dimension_number,1);

% This initializes a vector to store the size of the output arrays
array_size=zeros(1,array_dimension_number);

% This iterates through the different dimensions measuring the length of
% the vector in each dimension
for dimension_index=1:array_dimension_number;
    % This is the length of the current vector
    array_size(dimension_index)=length(vector_list{dimension_index});
end;

% This iterates through the different dimensions creating the gridded
% arrays
for dimension_index=1:array_dimension_number;
    % This is the vector of the current dimension
    current_vector=vector_list{dimension_index};
    
    % This is a permutation vector used to permute the current coordinate
    % vector into the correct dimension
    permutation_vector=1:array_dimension_number;
    % This sets the first value of the permutation vector to the current
    % dimension
    permutation_vector(1)=dimension_index;
    % This sets the dimension of the permutation vector to 1
    permutation_vector(dimension_index)=1;
    
    % This permutes the vector into the current dimension
    current_vector=permute(current_vector(:),permutation_vector);
    
    % This creates the current gridded array
    gridded_array{dimension_index}=bsxfun(@times,current_vector,ones(array_size));
    
end;



function linear_indices=sub2ind_vectorized(array_size,subscript_index_array);
% This function performs a similar operation to the builtin matlab function
% 'sub2ind' in that the function converts subscripted index values into
% linear index values.  However this function can be applied on vectors (of
% unknown length) which is beneficial for creating functions that can be
% called for an arbitrary number of dimensions.  The input argument
% 'array_size' is a 1 x N vector giving the size of the array into which
% the subscripts are indexed.  This array has N non-singleton dimensions.
% The input argument 'index_array' is a M x N array where M is the number
% of specific indices to be converted and N is the number of dimensions of
% the array.  The output argument 'linear_indices' is M x 1 in size and
% gives the linear index into the array for each of the M subscripted index
% sets.
%
% This function works by applying the following relationship
%
%  linear_index = ii + ( jj - 1 ) * ii_size + ( kk - 1 ) * ii_size * jj_size + . . .
%
% to calculate the linear indices.

% This is the vector of the products of the array dimensions
array_size_product=[1,cumprod(array_size(1:(end-1)))];

% This is the array giving the individual index components
subscript_addends=[subscript_index_array(:,1),subscript_index_array(:,2:end)-1];

% This takes the product of the array dimensions and the individual indices
subscript_product=bsxfun(@times,subscript_addends,array_size_product);

% This calculates the linear indices for the input subscripted indices by
% taking the sum of the individual products
linear_indices=sum(subscript_product,2);



function diagonal_matrix=create_diagonal_matrix(matrix_length,coefficient_kernal);
% This function creates a sparse diagonal matrix whose diagonal elements
% are given by the vector 'coefficient_kernal'.  The output sparse array is
% a 2D array that is 'matrix_length' x 'matrix_length' in size.

% This is 'half' of the length of the 'coefficient_kernal'
kernal_radius=(length(coefficient_kernal)-1)/2+1;

% This initializes the sparse diagonal matrix
diagonal_matrix=sparse([],[],[],matrix_length,matrix_length);

% This iterates through the coefficients adding them to the sparse array
for coefficient_index=1:length(coefficient_kernal);
    
    % This calculates the domain of the ii indices of the sparse array
    ii_vector_min=max([kernal_radius-coefficient_index+1,1]);
    ii_vector_max=matrix_length-max([coefficient_index-kernal_radius,0]);
    
    % This calculates the domain of the jj indices of the sparse array
    jj_vector_min=max([coefficient_index-kernal_radius,0])+1;
    jj_vector_max=matrix_length-max([kernal_radius-coefficient_index,0]);
    
    % This creates the vector of ii indices for the sparse array
    coefficient_ii_indices=(ii_vector_min:ii_vector_max)';
    % This creates the vector of jj indices for the sparse array
    coefficient_jj_indices=(jj_vector_min:jj_vector_max)';
    
    % This is the length of the current diagonal vector in the sparse array
    coefficient_vector_length=matrix_length-abs(kernal_radius-coefficient_index);
    
    % This creates the vector of coefficient values within the sparse array
    coefficient_values=coefficient_kernal(coefficient_index)*ones(coefficient_vector_length,1);
    
    % This adds the current set of coefficients to the initializes sparse
    % array
    diagonal_matrix=diagonal_matrix+sparse(coefficient_ii_indices,coefficient_jj_indices,coefficient_values,matrix_length,matrix_length);
    
end;



function x_round=precision_round(x,significant_decimal_number);
% This function rounds the array 'x' to the number of decimal points
% specified by 'significant_decimal_number'.  This is done in such a way 
% that the number of significant decimal points is the same regardless of
% the magnitude of the original array 'x', ie
%
%  precision_round(1e0*pi,4) -> 3.142
%  precision_round(1e2*pi,4) -> 314.2
%
% and so forth for other numbers.

% This is the scale factor to rescale 'x' so that it has an order of
% magnitude equal to 1
scale_factor=10.^-floor(log10(abs(x)));
% This rescales 'X' so that is has the specified precision
x_round=round(10^(significant_decimal_number-1)*scale_factor.*x)./(scale_factor*10^(significant_decimal_number-1));






	
	
	
	
