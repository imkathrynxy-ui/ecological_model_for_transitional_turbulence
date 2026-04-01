clearvars; clc; close all;
tic
advection_speed_range = [0.0165, 0.1125, 0.525];  % Food advection speed
n_ensemble = 2;
system_length = 3000;  % System length
system_width = 20;     % System width
threshold_density = 0.5; % Threshold to define turbulence
time_steps = 2e4;    % Number of time steps
plot_interval = 10;    % Interval for storing data
site_capacity = 5;     % Maximum food site capacity
V = 5;                 % System size (affects noise level)
t_nullcline = 1e4;     % Time to take snapshot for nullcline calculation

% save the location of front and back edge of slugs for velocity computation
front = zeros(time_steps/plot_interval,2,n_ensemble,length(advection_speed_range));
back = zeros(time_steps/plot_interval,2,n_ensemble,length(advection_speed_range));
nullcline = zeros(system_length,2,n_ensemble);
y_energy = zeros(system_length,7,10,n_ensemble);
for i = 1:length(advection_speed_range)

    food_adv_speed = advection_speed_range(i);  % Food advection speed

    % plot_data_x = zeros(time_steps/plot_interval,system_length);    % Storage predator density data for plot
    % plot_data_y = zeros(time_steps/plot_interval,system_length);    % Storage prey density data for plot
    % plot_data_g = zeros(time_steps/plot_interval,system_length);    % Storage nutrient density data for plot
    %
    for e = 1:n_ensemble

        D = 0.125;           % diffusion rate  (paper: D_A = D_B = 0.125)
        d1 = 0.02;           % death rate of predator  (paper: d_A = 0.02)
        d2 = 0.02;           % death rate of prey      (paper: d_B = 0.02)
        p = 0.05;            % predation rate          (paper: p = 0.05)
        m = 0.0002;          % mutation rate           (paper: m = 0.0002)
        b = 0.05;            % birth rate of prey      (paper: b = 0.05)
        c = 0.04;            % competition rate        (paper: c_A = c_B = 0.04)
        growth_rate = 20/9 * food_adv_speed^2; % growth rate of grass (paper: g = 20U^2/9)

        % Rescaling all parameters by food advection speed
        D = rescale_by_speed(D,food_adv_speed);
        d1 = rescale_by_speed(d1,food_adv_speed);
        d2 = rescale_by_speed(d2,food_adv_speed);
        p = rescale_by_speed(p,food_adv_speed);
        m = rescale_by_speed(m,food_adv_speed);
        b = rescale_by_speed(b,food_adv_speed);
        c = rescale_by_speed(c,food_adv_speed);
        growth_rate = rescale_by_speed(growth_rate,food_adv_speed);

        % Storage data for space-time plots
        plot_data_x = zeros(time_steps/plot_interval, system_length);
        plot_data_y = zeros(time_steps/plot_interval, system_length);
        plot_data_g = zeros(time_steps/plot_interval, system_length);

        % Initialize population grids
        x = zeros(system_width, system_length); % Predators
        y = zeros(system_width, system_length); % Prey
        g = site_capacity * ones(system_width, system_length); % Food

        % Initial conditions
        mid = system_length - 1500;
        x(:, mid-15:mid+15) = rand(system_width, 31) > 2/5;
        y(:, mid-15:mid+15) = rand(system_width, 31) > 2/5;
        g(:, mid+15:end) = 0;

        % Time evolution
        for t = 1:time_steps

            
            y_loss_p = zeros(system_width,system_length);
            y_loss_d = zeros(system_width,system_length);
            y_loss_c = zeros(system_width,system_length);
            y_gain_b = zeros(system_width,system_length);
            y_loss_m = zeros(system_width,system_length);
            y_diffuse = zeros(system_width,system_length);
            x_diffuse = zeros(system_width,system_length);
            y_last_step = y;  % snapshot of y at the start of this step

            g(:,1) = site_capacity; % Refill food on left boundary

            for k = 1:system_width * system_length
                jx = randi(system_width);
                jy = randi(system_length);
                s = rand();

                if s < 1/10
                    % Diffusion
                    [x, x_diffuse] = diffuse(x, jx, jy, D, system_width, system_length, x_diffuse);
                elseif s < 2/10
                    [y, y_diffuse] = diffuse(y, jx, jy, D, system_width, system_length, y_diffuse);
                elseif s < 3/10
                    [y, g, y_gain_b] = reproduce(y, g, jx, jy, b, V, system_width, system_length, y_gain_b);
                elseif s < 4/10
                    [x, y, y_loss_p] = predation(x, y, jx, jy, p, V, system_width, system_length, y_loss_p);
                elseif s < 5/10
                    if rand() < 1 - exp(-d1 * x(jx, jy))
                        x(jx, jy) = max(0, x(jx, jy) - 1);
                    end
                elseif s < 6/10
                    if rand() < 1 - exp(-d2 * y(jx, jy))
                        y(jx, jy) = max(0, y(jx, jy) - 1);
                        y_loss_d(jx,jy)=y_loss_d(jx,jy)-1;
                    end
                elseif s < 7/10
                    [x, y, y_loss_m] = mutate(x, y, jx, jy, m, y_loss_m);
                elseif s < 8/10
                    if rand() < 1 - exp(-c * x(jx, jy) * (x(jx, jy) - 1) / V)
                        x(jx, jy) = max(0, x(jx, jy) - 1);
                    end
                elseif s < 9/10
                    if rand() < 1 - exp(-c * y(jx, jy) * (y(jx, jy) - 1) / V)
                        y(jx, jy) = max(0, y(jx, jy) - 1);
                        y_loss_c(jx,jy)=y_loss_c(jx,jy)-1;
                    end
                else
                    if rand() < 1 - exp(-growth_rate * V)
                        g(jx, jy) = min(site_capacity, g(jx, jy) + 1);
                    end
                end
            end

            % Food advection
            g = circshift(g, 1, 2);

            % Store data for plotting
            if mod(t, plot_interval) == 0
                plot_data_x(t / plot_interval, :) = mean(x);  %average over the width of the pipe
                plot_data_y(t / plot_interval, :) = mean(y);
                plot_data_g(t / plot_interval, :) = mean(g);

                pos = find(mean(y)>=threshold_density);    % find out the position of turbulent sites

                if length(pos)>1
                    front(t/plot_interval,:,e,i) = [t/food_adv_speed,pos(1)];
                    back(t/plot_interval,:,e,i) = [t/food_adv_speed,pos(end)];
                end

            end

            if t == t_nullcline
                nullcline(:,1,e) = mean(g)';
                nullcline(:,2,e) = mean(y)';
            end

            
            if t == time_steps
                y_energy(:,1,time_steps-t+1,e) = mean(y_gain_b)';   % 1: gain from birth
                y_energy(:,2,time_steps-t+1,e) = mean(y_loss_p)'; % 2: loss from predation
                y_energy(:,3,time_steps-t+1,e) = mean(y_loss_c)'; % 3: loss from competition
                y_energy(:,4,time_steps-t+1,e) = mean(y_loss_d)';  % 4: loss from death
                y_energy(:,5,time_steps-t+1,e) = mean(y_loss_m)';  % 5: loss from mutation
                y_energy(:,6,time_steps-t+1,e) = mean(y_diffuse)';  %6: gain and loss from diffusion
                y_energy(:,7,time_steps-t+1,e) = mean(y_last_step-y)';  %net gain in y
            end

        end
        y_energy_save = sum(y_energy,3);
   end
end
[front_velocity_mean,front_velocity_std,back_velocity_mean,back_velocity_std] = calculate_front_velocity(front, back, advection_speed_range, plot_interval);
mean_nullcline = plot_nullcline(nullcline);
plot_y_energy(y_energy);

toc

% Function definitions
function y = rescale_by_speed(x,speed)
    y = x / speed;
end

function [grid, y_diffuse] = diffuse(grid, jx, jy, D, width, length, y_diffuse)
    if rand() < 1 - exp(-D * grid(jx, jy))
        dir = randi(4);
        [jx_new, jy_new] = get_neighbor(jx, jy, dir, width, length);
        grid(jx_new, jy_new) = grid(jx_new, jy_new) + 1;
        grid(jx, jy) = grid(jx, jy) - 1;

        y_diffuse(jx_new,jy_new) = y_diffuse(jx_new,jy_new)+1;
        y_diffuse(jx,jy) = y_diffuse(jx,jy)-1;
    end
end

function [grid, food, y_gain_b] = reproduce(grid, food, jx, jy, b, V, width, length, y_gain_b)
    dir = randi(4);
    [jx_new, jy_new] = get_neighbor(jx, jy, dir, width, length);
    if rand() < 1 - exp(-b * grid(jx, jy) * food(jx_new, jy_new) / V)
        grid(jx_new, jy_new) = grid(jx_new, jy_new) + 1;
        food(jx_new, jy_new) = max(0, food(jx_new, jy_new) - 1);

        y_gain_b(jx_new,jy_new)=y_gain_b(jx_new,jy_new)+1;
    end
end

function [x, y, y_loss_p] = predation(x, y, jx, jy, p, V, width, length, y_loss_p)
    dir = randi(4);
    [jx_new, jy_new] = get_neighbor(jx, jy, dir, width, length);
    if rand() < 1 - exp(-p * x(jx, jy) * y(jx_new, jy_new) / V)
        x(jx_new, jy_new) = x(jx_new, jy_new) + 1;
        y(jx_new, jy_new) = max(0, y(jx_new, jy_new) - 1);
        y_loss_p(jx_new, jy_new) = y_loss_p(jx_new, jy_new) - 1;
    end
end

function [x, y, y_loss_m] = mutate(x, y, jx, jy, m, y_loss_m)
    if rand() < 1 - exp(-m * y(jx, jy))
        x(jx, jy) = x(jx, jy) + 1;
        y(jx, jy) = max(0, y(jx, jy) - 1);
        y_loss_m(jx,jy)=y_loss_m(jx,jy)-1;
    end
end

function [jx_new, jy_new] = get_neighbor(jx, jy, dir, width, length)
    switch dir  %periodic boundary conditions for all four boundaries
        case 1, jx_new = mod(jx - 2, width) + 1; jy_new = jy;   % Move up
        case 2, jx_new = mod(jx, width) + 1; jy_new = jy;       % Move down
        case 3, jy_new = mod(jy - 2, length) + 1; jx_new = jx;  % Move left
        case 4, jy_new = mod(jy, length) + 1; jx_new = jx;      % Move right
    end
end

function plot_space_time_data(system_length, time_steps, plot_interval, u, plot_data_x, plot_data_y, plot_data_g)
    [X,Y] = meshgrid(1:system_length,1:time_steps/plot_interval);
    Y = Y*plot_interval/u;
    figure
    subplot(3,1,1)
    contourf(X,Y,plot_data_x(:,:,1),'linestyle','none')
    xlabel('space')
    ylabel('time')
    title('predator')
    subplot(3,1,2)
    contourf(X,Y,plot_data_y(:,:,1),'linestyle','none')
    hold on
    xlabel('space')
    ylabel('time')
    title('prey')
    subplot(3,1,3)
    contourf(X,Y,plot_data_g(:,:,1),'linestyle','none')
    xlabel('space')
    ylabel('time')
    title('nutrient')
end

function [front_velocity_mean,front_velocity_std,back_velocity_mean,back_velocity_std] = calculate_front_velocity(front, back, advection_speed_range, plot_interval)
    [num_time_steps, ~, n_ensemble, num_speeds] = size(front);

    % Initialize storage for velocity statistics
    front_velocity_mean = zeros(num_speeds, 1);
    front_velocity_std = zeros(num_speeds, 1);
    back_velocity_mean = zeros(num_speeds, 1);
    back_velocity_std = zeros(num_speeds, 1);

    % Loop over advection speeds
    for i = 1:num_speeds
        front_velocities = zeros(n_ensemble,1);
        back_velocities = zeros(n_ensemble,1);

        for e = 1:n_ensemble
            time_diff = diff(squeeze(front(:,1,e,i))) * plot_interval; % Compute time differences
            space_diff_front = diff(squeeze(front(:,2,e,i))); % Compute space differences
            space_diff_back = diff(squeeze(back(:,2,e,i)));

            valid_idx = time_diff > 0; % Ensure valid division

            % Compute velocities
            velocity_front = space_diff_front(valid_idx) ./ time_diff(valid_idx);
            velocity_back = space_diff_back(valid_idx) ./ time_diff(valid_idx);

            front_velocities(e) = mean(velocity_front);
            back_velocities(e) = mean(velocity_back);
        end

        % Compute mean and std
        front_velocity_mean(i) = mean(front_velocities);
        front_velocity_std(i) = std(front_velocities);
        back_velocity_mean(i) = mean(back_velocities);
        back_velocity_std(i) = std(back_velocities);
    end

    % Plot results
    figure;
    errorbar(advection_speed_range, front_velocity_mean, front_velocity_std, 'o', 'LineWidth', 1.5);
    hold on;
    errorbar(advection_speed_range, back_velocity_mean, back_velocity_std, 's', 'LineWidth', 1.5);
    hold off;
    xlabel('Advection Velocity');
    ylabel('Frint speed');
    legend('Upstream Front', 'Downstream Front');
    title('Velocity of Turbulence Front and Back Edge');
end

function mean_nullcline = plot_nullcline(nullcline)
    % plot_nullcline: Average over ensembles and plot the nullcline.
    %
    % Input:
    %   nullcline: a matrix of size (system_length x 2 x n_ensemble)
    %              The first column is the nutrient (g) profile,
    %              and the second column is the prey (y) profile.
    %
    % The function computes the mean nullcline over the ensembles and
    % plots both nutrient and prey density as functions of spatial position.

    % Determine the spatial size and number of ensembles
    [system_length, numVars, n_ensemble] = size(nullcline);

    % Check that we have exactly two variables (g and y)
    if numVars ~= 2
        error('The second dimension of nullcline must be 2 (nutrient and prey)');
    end

    % Average over ensembles
    mean_nullcline = mean(nullcline, 3);  % size: system_length x 2

    % Create spatial coordinate (assuming unit spacing)
    x_axis = 1:system_length;

    % Plotting
    figure;
    subplot(2,1,1)
    plot(x_axis, mean_nullcline(:,1), 'g-', 'LineWidth', 2);
    xlabel('Spatial Position');
    ylabel('Nutrient Density');
    title('Averaged Nullcline Profiles');
    grid on;

    subplot(2,1,2)
    plot(x_axis, mean_nullcline(:,2), 'b-', 'LineWidth', 2);
    xlabel('Spatial Position');
    ylabel('Prey Density');
    title('Averaged Nullcline Profiles');
    grid on;

    figure;
    plot(mean_nullcline(:,1), mean_nullcline(:,2), '-', 'LineWidth', 2);
    xlabel('Nurient density');
    ylabel('Prey density');
end

function plot_y_energy(y_energy)
    % plot_y_energy: Plot the spatial profiles of different energy contributions
    % for the prey dynamics.
    %
    % Input:
    %   y_energy: a 4D array of size (system_length x 7 x T x n_ensemble)
    %             where the second dimension corresponds to:
    %               1. Gain from birth
    %               2. Loss from predation
    %               3. Loss from competition
    %               4. Loss from death
    %               5. Loss from mutation
    %               6. Contribution from diffusion (gain/loss)
    %               7. Net gain in y
    %
    % The third dimension is (optionally) time or another index.
    % The function first sums (or averages) over the third dimension and over ensembles,
    % then plots the resulting profiles for each contribution.

    % Get the size of the data
    [system_length, numContributions, T, n_ensemble] = size(y_energy);

    % For this plot, we sum over the third dimension (time or other index)
    % and average over ensembles.
    y_energy_sum = sum(y_energy, 3);  % now size: (system_length x 7 x n_ensemble)
    y_energy_avg = mean(y_energy_sum, 3);  % now size: (system_length x 7)

    % Create spatial coordinate
    x_axis = 1:system_length;

    % Define labels for the contributions
    energy_labels = {'Birth Gain', 'Predation Loss', 'Competition Loss', ...
                     'Death Loss', 'Mutation Loss', 'Diffusion', 'Net Gain'};

    % Plot each contribution in a subplot or as multiple curves in one figure.
    % Here, we plot all contributions in one figure for comparison.
    figure;
    hold on;
    colors = lines(numContributions);
    for idx = 1:numContributions
        plot(x_axis, y_energy_avg(:, idx), 'Color', colors(idx,:), 'LineWidth', 2);
        hold on;
    end
    xlabel('Spatial Position');
    ylabel('Energy Contribution (per time step)');
    title('Spatial Profiles of Energy Gain and Loss for Prey (y)');
    legend(energy_labels, 'Location', 'Best');
    grid on;

    figure;
    subplot(3,1,1)
    plot(x_axis, y_energy_avg(:,1), 'LineWidth', 2)
    xlabel('Spatial Position')
    ylabel('Rate of prey density gain')
    subplot(3,1,2)
    plot(x_axis, sum(y_energy_avg(:,2:5),2), 'LineWidth', 2)
    xlabel('Spatial Position')
    ylabel('Rate of prey density loss')
    subplot(3,1,3)
    plot(x_axis, y_energy_avg(:,7), 'LineWidth', 2)
    xlabel('Spatial Position')
    ylabel('Rate of prey density net change')
end
