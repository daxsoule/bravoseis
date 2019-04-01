function createLines
% function to create a bunch of shooting lines

% Load the OBSs
% load ../final/GE-OBS.xy	
% load ../final/US-OBS.xy	

% Distance scaling
lat2km = 111.3;
lon2km = lat2km * cosd(62.5);

%% 3D Profiles & Tomography
mid3d = [-59.8408 -62.8541];     % Center point of 3D survey
mid3d = [-59.8408 -62.8541]     % Center point of 3D survey
theta3d = 50;                   % Azimuth of shot lines
len3d = 15;                      % Length of shot lines   
addLine = 0;                    % Distance to add at either end of 3d lines
spacing3d = 1.8;                 % Spacing of shot lines
n3d = 7; %41                       % Number of 3d lines
spacing3dTie = 3.75; %3.75;             % Spacing of the tie lines
n3dTie = 3;                      % Number of tie lines
%lenTomo = 15;                    % Length of tomography lines
%nTomo = 30;                      % Number of tomography lines
%taperTomo = [3 2.5 2 1.5 1 0.5]; % Tapering outer tomography lines
speed = 4.0*1.85;                % shooting speed in km/hr 
radiusTurn = 1.5;                % turning radius
radiusTurnTomo = 0.25;           % Turning radius for tomography

% Lines
x0 = (-floor(n3d/2):floor(n3d/2))*spacing3d;
x1 = x0;
y0 = -ones(1,n3d)*(len3d/2+addLine);
y1 = -y0;
x0Tie = -ones(1,n3dTie)*spacing3d*floor(n3d/2);
x1Tie = -x0Tie;
y0Tie = -floor(n3dTie/2)*spacing3dTie:spacing3dTie:floor(n3dTie/2)*spacing3dTie;
y1Tie = y0Tie;
%x0Tomo = (x0(1:end-1)+x0(2:end))/2;
%x1Tomo = (x1(1:end-1)+x1(2:end))/2;
%i = 1+(n3d-nTomo-1)/2:(n3d-1)-(n3d-nTomo-1)/2;
% x0Tomo = x0Tomo(i);
% x1Tomo = x1Tomo(i);
% y0Tomo = -ones(1,nTomo)*lenTomo/2;
% taper = zeros(size(y0Tomo));
% taper(1:length(taperTomo)) = taperTomo;
% taper(end:-1:end-length(taperTomo)+1) = taperTomo;
% y0Tomo = y0Tomo + taper;
% y1Tomo = -y0Tomo;

% Time calculations
t3d = (sum(abs(y1-y0))+pi*radiusTurn*(n3d-1))/speed
t3dtie =  (sum(abs(x1Tie-x0Tie))+pi*spacing3dTie/2*(length(x1Tie)-1))/speed
%tTomo = (sum(abs(y1Tomo-y0Tomo))+pi*radiusTurnTomo*(nTomo-1))/speed

%Rotate
x0p = (cosd(theta3d)*x0 + sind(theta3d)*y0)/lon2km+mid3d(1);
y0p = (-sind(theta3d)*x0 + cosd(theta3d)*y0)/lat2km+mid3d(2);
x1p = (cosd(theta3d)*x1 + sind(theta3d)*y1)/lon2km+mid3d(1);
y1p = (-sind(theta3d)*x1 + cosd(theta3d)*y1)/lat2km+mid3d(2);
% x0Tomop = (cosd(theta3d)*x0Tomo + sind(theta3d)*y0Tomo)/lon2km+mid3d(1);
% y0Tomop = (-sind(theta3d)*x0Tomo + cosd(theta3d)*y0Tomo)/lat2km+mid3d(2);
% x1Tomop = (cosd(theta3d)*x1Tomo + sind(theta3d)*y1Tomo)/lon2km+mid3d(1);
% y1Tomop = (-sind(theta3d)*x1Tomo + cosd(theta3d)*y1Tomo)/lat2km+mid3d(2);
x0Tiep = (cosd(theta3d)*x0Tie + sind(theta3d)*y0Tie)/lon2km+mid3d(1);
y0Tiep = (-sind(theta3d)*x0Tie + cosd(theta3d)*y0Tie)/lat2km+mid3d(2);
x1Tiep = (cosd(theta3d)*x1Tie + sind(theta3d)*y1Tie)/lon2km+mid3d(1);
y1Tiep = (-sind(theta3d)*x1Tie + cosd(theta3d)*y1Tie)/lat2km+mid3d(2);

% PLot
close all
%plot(GE_OBS(:,1),GE_OBS(:,2),'sk','markerfacecolor','k')

%plot(US_OBS(:,1),US_OBS(:,2),'sk','markerfacecolor','k')
plot([x0p; x1p],[y0p; y1p],'-k')
hold on
plot([x0Tiep; x1Tiep],[y0Tiep; y1Tiep],'-k')
%plot([x0Tomop; x1Tomop],[y0Tomop; y1Tomop],'-r')

fid = fopen('y_bathy_a.txt','wt');
for i=1:length(x0p)
  fprintf(fid,'%f %f\n',x0p(i),y0p(i));
  fprintf(fid,'%f %f\n',x1p(i),y1p(i));
  if i<length(x0p)
    fprintf(fid,'>\n');
  end
end
fid = fopen('x_bathy_a.txt','wt');
for i=1:length(x0Tiep)
  fprintf(fid,'%f %f\n',x0Tiep(i),y0Tiep(i));
  fprintf(fid,'%f %f\n',x1Tiep(i),y1Tiep(i));
  if i<length(x0Tiep)
    fprintf(fid,'>\n');
  end
end
fclose(fid);
% fid = fopen('tomo_lines.txt','wt');
% for i=1:length(x0Tomop)
%   fprintf(fid,'%f %f\n',x0Tomop(i),y0Tomop(i));
%   fprintf(fid,'%f %f\n',x1Tomop(i),y1Tomop(i));
%   if i<length(x0Tomop)
%     fprintf(fid,'>\n');
%   end
% end
% fclose(fid);


%% 2D Profiles
mid2d = [-59 -62.6667];
len2d = 20;
theta2dLeft = -31.5;
theta2dRight = -31.5;
xoffset2d = -25;
yoffset2dLeft = 4;
yoffset2dRight = 1.5;
spacing2d = 5;
nleft = round(60/spacing2d);
nright = round(24/spacing2d);

x0 = (-nleft:nright)*spacing2d;
y0 = -ones(1,nleft+nright+1)*(len2d/2+addLine);
x1 = x0;
y1 = -y0;

i = x0<xoffset2d;
y0(i) = y0(i)+yoffset2dLeft;
y1(i) = y1(i)+yoffset2dLeft;
x02dp(i) = (cosd(theta2dLeft)*x0(i) + sind(theta2dLeft)*y0(i))/lon2km+mid2d(1);
y02dp(i) = (-sind(theta2dLeft)*x0(i) + cosd(theta2dLeft)*y0(i))/lat2km+mid2d(2);
x12dp(i) = (cosd(theta2dLeft)*x1(i) + sind(theta2dLeft)*y1(i))/lon2km+mid2d(1);
y12dp(i) = (-sind(theta2dLeft)*x1(i) + cosd(theta2dLeft)*y1(i))/lat2km+mid2d(2);
i = x0>=xoffset2d;
y0(i) = y0(i)+yoffset2dRight;
y1(i) = y1(i)+yoffset2dRight;
x02dp(i) = (cosd(theta2dRight)*x0(i) + sind(theta2dRight)*y0(i))/lon2km+mid2d(1);
y02dp(i) = (-sind(theta2dRight)*x0(i) + cosd(theta2dRight)*y0(i))/lat2km+mid2d(2);
x12dp(i) = (cosd(theta2dRight)*x1(i) + sind(theta2dRight)*y1(i))/lon2km+mid2d(1);
y12dp(i) = (-sind(theta2dRight)*x1(i) + cosd(theta2dRight)*y1(i))/lat2km+mid2d(2);


% Time
t2d = (sum(abs(y1-y0))+(pi*radiusTurn+spacing2d-2*radiusTurn)*(length(y0)-1))/speed


plot([x02dp; x12dp],[y02dp; y12dp],'-g')

fid = fopen('2D_lines.txt','wt');
for i=1:length(x02dp)
  fprintf(fid,'%f %f\n',x02dp(i),y02dp(i));
  fprintf(fid,'%f %f\n',x12dp(i),y12dp(i));
  if i<length(x02dp)
    fprintf(fid,'>\n');
  end
end
fclose(fid);

%% 2d Along Axis line
data = [-60.07184	-62.90529	-1002.8	0	0.00
-59.87934	-62.85347	-589.1	11.33	11.33
-59.79500	-62.82920	-626.3	16.39	5.06
-59.71616	-62.81831	-983.7	20.58	4.18
-59.63366	-62.79904	-1103.3	25.28	4.71
-59.43198	-62.77641	-1434.9	35.84	10.56
-59.27064	-62.73024	-1372.4	45.53	9.69
-59.21381	-62.71428	-1369.5	48.93	3.40
-59.05797	-62.66970	-1242.4	58.29	9.37
-58.98280	-62.64275	-1076.4	63.16	4.87
-58.87279	-62.61999	-1249.1	69.33	6.17
-58.74445	-62.57611	-1415.4	77.51	8.18
-58.65278	-62.53640	-1646.4	83.96	6.45
-58.595399	-62.497371	    	-1334.9	94.08	10.12
-58.284601	-62.389229	-1488.4	107.30	13.23];

lon = data(:,1);
lat = data(:,2);
lon(end-1,1) = x0Tiep(ceil(n3dTie/2));
lon(end,1) = x1Tiep(ceil(n3dTie/2));
lat(end-1,1) = y0Tiep(ceil(n3dTie/2));
lat(end,1) = y1Tiep(ceil(n3dTie/2));

plot(lon,lat,'-g')
xlim([-60.2 -58.1])
ylim([-63 -62.3])

% Time
tAlong = sum(sqrt((diff(lon)*111.19*cosd(62.5)).^2+(diff(lat)*111.19).^2))/speed;

fid = fopen('alongAxis_line.txt','wt');
for i=1:length(lon)
  fprintf(fid,'%f %f\n',lon(i),lat(i));
end
fclose(fid);

%% Write Excel File
fid = fopen('ShootingPlan.csv','wt+');
fprintf(fid,'\nAlong Axis Line 101\n');
fprintf(fid,'Way Point,Lat,Lon,Lat Deg,Lat Min,Lon Deg,Lon Min\n');
data = [(1:length(lat))' lat lon floor(-lat) 60*(-lat-floor(-lat)) floor(-lon) 60*(-lon-floor(-lon))];
fprintf(fid,'%i, %f, %f, %i, %f, %i, %f\n',data');
write_lines(fid,'2-D Lines',200,y02dp,y12dp,x02dp,x12dp)
write_lines(fid,'3-D Reflection Lines',300,y0p,y1p,x0p,x1p)
write_lines(fid,'3-D Tie Lines',400,y0Tiep,y1Tiep,x0Tiep,x1Tiep)
write_lines(fid,'Tomography Lines',500,y0Tomop,y1Tomop,x0Tomop,x1Tomop)
fclose(fid);

end

%%
function [deg,min] = deg2degmin(decdeg)
% Converts decimal degrees to degrees and minutes (unsigned)

% Decimal places for minutes
dp = 4;

decdeg = abs(decdeg);
deg = floor(decdeg);
min = 60*(decdeg-deg);
min = round(min*10^dp)/10^dp;

end

%%
function write_lines(fid,titl,linenumber,y0,y1,x0,x1)

fprintf(fid,['\n\n' titl '\n']);
fprintf(fid,'Line,,Waypoint,Lat,Lon,Lat Deg,Lat Min,Lon Deg,Lon Min,,Waypoint,Lat,Lon,Lat Deg,Lat Min,Lon Deg,Lon Min\n');
[y0d,y0m] = deg2degmin(y0);
[y1d,y1m] = deg2degmin(y1);
[x0d,x0m] = deg2degmin(x0);
[x1d,x1m] = deg2degmin(x1);
for i = 1:length(x0)
  w0 = [int2str(linenumber+i) 'A'];
  w1 = [int2str(linenumber+i) 'B'];
  fprintf(fid,'%i,, %s,%f, %f, %i, %f, %i, %f, ,%s ,%f, %f, %i, %f, %i, %f\n', linenumber+i, ...
                 w0,y0(i),x0(i),y0d(i),y0m(i),x0d(i),x0m(i), ...
                 w1,y1(i),x1(i),y1d(i),y1m(i),x1d(i),x1m(i));
end

end
