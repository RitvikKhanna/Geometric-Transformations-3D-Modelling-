clear all;
close all;

% Generate random spherical coordinates(radial, azimuthal, polar)
N = 5000;
radial = 8*(rand(N,1)); %
azimuthal = 2*pi*rand(N,1); %azimuthal
polar = pi*rand(N,1); %radial

% Convert to cartesian
x = radial.*cos(azimuthal).*sin(polar);
y = radial.*sin(azimuthal).*sin(polar);
z = radial.*cos(polar);
% Convert to Homogeneous Coordinate
sphere = [x y z ones(5000,1)];
% I do this ones(5000,1) to convert them to homogenous by adding one extra
% point whose value is 1, for every point

% 3D Translation, on Original Sphere
% Random Translation by (2 3 5)
% For truly random values we can use 3 variables
translation1 = 35.*(rand(1));
% similarly for the other 2 transations translation2 & translation3
translation2 = 35.*(rand(1));
translation3 = 35.*(rand(1));
% Then substitue values translation1 instead of 75 below, translation2
% instead of 25, and translation3 instead of 50
M1 = [1 0 0 translation1; 0 1 0 translation2; 0 0 1 translation3; 0 0 0 1];
size(M1)

% sphere_M1 = ;
% I do this since we have 5000 points so I create a 5000x4 Ones matrix
sphere_M1 = ones(5000,4);

% Now I use a for loop to translate each and every point of those 5000
% points of the sphere and then put them in sphere_M1
% For the row starting at row 1 and incrementing after all column iterations are done
for row=1:5000
    % Do it for all the 4 columns
    result = M1*[sphere(row, 1); sphere(row, 2); sphere(row, 3); sphere(row, 4)];
    % Now for all the 4 rows store it the translated point in sphere_M1
    for col=1:4
        sphere_M1(row, col) = result(col);
    end
end


% 3D Translation with Non-uniform scaling, on Original Sphere
% For truly random values we can use 3 variables
scale1 = (2.5).*(rand(1));
% similarly for the other 2 transations scale2 & scale3
scale2 = (2.5).*(rand(1));
scale3 = (2.5).*(rand(1));
% Then substitue values scale1 instead of 0.5 below, scale2 instead of 1.5, 
% and scale3 instead of 0.8
M2 = [scale1 0 0 0; 0 scale2 0 0; 0 0 scale3 0; 0 0 0 1];
size(M2)

% Using the same logic as above and using a for loop to iterate through all
% the 5000 points i.e. 5000 rows and 4 columns.
% We do for every point M2*M1.
% sphere_M2 = ;
sphere_M2 = ones(5000,4);

for row=1:5000
    % M2*M1*Point because we have translation followed by scaling
    result = M2*M1*[sphere(row, 1); sphere(row, 2); sphere(row, 3); sphere(row, 4)];
    for col=1:4
        sphere_M2(row, col) = result(col);
    end
end

% 3D 11Rotations on Sphere M2
% Using a random angle as 55 degrees.
% We can generate a random angle theta
theta = 180.*rand(1);
% theta1 rotation along x axis
theta1 = [1 0 0 0; 0 cosd(theta) -sind(theta) 0; 0 sind(theta) cos(theta) 0; 0 0 0 1];
% theta2 rotation along y axis
theta2 = [cosd(theta) 0 sind(theta) 0; 0 1 0 0; -sind(theta) 0 cosd(theta) 0; 0 0 0 1];
%theta3 = rotation along z axis
theta3 = [cosd(theta) -sind(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];

%We declare the final roation among all axis's as 
rotation = theta3*theta2*theta1;

% We use the same logic of 5000 points i.e. 5000x4 matrix for each point
% roation
sphere_R = ones(5000,4);

% We use a for row 1-5000 and for every column 1-4 rotation of all the
% points using our final roation defined as 
% rotation = theta3*theta2*theta1;

for row=1:5000
    result = rotation*[sphere_M2(row, 1); sphere_M2(row, 2); sphere_M2(row, 3); sphere_M2(row, 4)];
    for col=1:4
        sphere_R(row, col) = result(col);
    end
end


% Composite Transformation
sphere_M = ones(5000,4);

% Let our composite matrix be
composite = rotation*M2*M1;

for row=1:5000
    result = composite*[sphere(row, 1); sphere(row, 2); sphere(row, 3); sphere(row, 4)];
    for col=1:4
        sphere_M(row, col) = result(col);
    end
end

% Check the norm after applying two equivalent transformations
fprintf('Norm: %i\n', norm(sphere_R - sphere_M))

% plot
figure(1)
hold on;
scatter3(x, y, z, 2, 'filled', 'r')
scatter3(sphere_M1(:,1), sphere_M1(:,2), sphere_M1(:,3), 2, 'filled', 'b')
scatter3(sphere_M2(:,1), sphere_M2(:,2), sphere_M2(:,3), 2, 'filled', 'y')
scatter3(sphere_R(:,1), sphere_R(:,2), sphere_R(:,3), 2, 'filled', 'g')
scatter3(sphere_M(:,1), sphere_M(:,2), sphere_M(:,3), 2, 'filled', 'p')
axis vis3d equal
grid on;
xlabel X, ylabel Y, zlabel Z
hold off;