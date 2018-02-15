clear all
close all
clc


%%%%%%%%%%%%%%%%%%
%%% Parameters %%%
%%%%%%%%%%%%%%%%%%

di_const1 = 1;
di_const2 = 12;
% di_const2 = 2;

hmax = 2 * pi * 10 / 20;

lambda = 700;

k0 = 2*pi/lambda;

%% Cylinder specifics

n_cyl = 7;
r_cyl = 10;
cyl_period = 50;
n_points_cyl = 60;

%% Area specifics

ul_spacing = 1500;
area_period = 30;
area_height = n_cyl * 300;
tot_height = ul_spacing + area_height;
n_points_ud = 50;
n_points_lr = floor((ul_spacing + area_height)/10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Structure construction %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Centres of cylinders.

if mod(n_cyl,2) == 0				% Checks if n_cyl is an even number.
	
	cyl_cent_x = 0 * ones(1,n_cyl);
	cyl_cent_y = zeros(1,n_cyl);
	
	for i = 1:n_cyl/2
		
		cyl_cent_y(2*i-1) = i * cyl_period - cyl_period/2;
		cyl_cent_y(2*i) = - i * cyl_period + cyl_period/2;
		
	end
	
else
	
	cyl_cent_x = 0 * ones(1,n_cyl);
	cyl_cent_y = zeros(1,n_cyl);
	
	for i = 1:(n_cyl-1)/2
		cyl_cent_y(2*i) = i * cyl_period;
		cyl_cent_y((2*i)+1) = - i * cyl_period;
	end
	
end

%{
%% Points on cylinders' periphery.

cyl_peri_x = zeros(1,n_points_cyl);
cyl_peri_y = zeros(1,n_points_cyl);
cyl_cell = cell(2,n_cyl);

for j = 1:n_cyl
	
	for i = 1:n_points_cyl
		
		cyl_peri_x(i+n_points_cyl*(j-1)) = cyl_cent_x(j) + cos(2*pi/n_points_cyl*i) * r_cyl;
		cyl_peri_y(i+n_points_cyl*(j-1)) = cyl_cent_y(j) + sin(2*pi/n_points_cyl*i) * r_cyl;
		
	end
	
	cyl_cell{1,j} = cyl_peri_x(((j-1)*n_points_cyl)+1:j*n_points_cyl);
	cyl_cell{2,j} = cyl_peri_y(((j-1)*n_points_cyl)+1:j*n_points_cyl);
	
end

%% Points on Area Boundary

area_points_x = zeros(1,2 * n_points_lr + 2 * n_points_ud);
area_points_y = zeros(1,2 * n_points_lr + 2 * n_points_ud);

counter = 0;

for l = 1:length(area_points_x)
	
	counter = counter + 1;
	
	if l <= n_points_ud
		
		key = 1;
    
        area_points_x(l) = -area_period/2 + ((area_period / n_points_ud) * counter);

        area_points_y(l) = -(area_height/2 + ul_spacing);
		
	elseif key == 1
		
		counter = 1;
		
	end
	
	if n_points_ud < l && l <= n_points_ud + n_points_lr
		
		key = 2;
		
		area_points_x(l) = area_period/2;
		area_points_y(l) = -(area_height/2 + ul_spacing) + (((2 * ul_spacing + area_height) / n_points_lr) * counter);
		
	elseif key == 2
		
		counter = 1;
    
	end
	
	if n_points_ud + n_points_lr < l && l <= 2 * n_points_ud + n_points_lr
		
		key = 3;
		
		area_points_x(l) = area_period/2 - ((area_period / n_points_ud) * counter);
		area_points_y(l) = area_height/2 + ul_spacing;
		
	elseif key == 3
		
		counter = 1;
		
	end
	
	if 2 * n_points_ud + n_points_lr < l && l <= 2 * n_points_ud + 2 * n_points_lr
		
		key = 4;
		
		area_points_x(l) = -area_period/2;
		area_points_y(l) = area_height/2 + ul_spacing - ((2 * ul_spacing + area_height)/n_points_lr * counter);
		
	elseif key == 4
		
		counter = 1;
		
	end
	
end

%% Create Geometry Matrix

geom = zeros(length(area_points_x)+n_points_cyl,7);

for i = 1:length(area_points_x)-1
	
	geom(i,:) = [2 area_points_x(i) area_points_x(i+1) area_points_y(i) area_points_y(i+1) 2 0];
	
end

geom(length(area_points_x),:) = [2 area_points_x(length(area_points_x)) area_points_x(1) area_points_y(length(area_points_y)) area_points_y(1) 2 0];

var_len = length(area_points_x);

for j = 1:n_cyl
	
	for i = 1:n_points_cyl-1
		
		geom(i+var_len,:) = [2 cyl_cell{1,j}(i) cyl_cell{1,j}(i+1) cyl_cell{2,j}(i) cyl_cell{2,j}(i+1) 1 2];%[2 cyl_peri_x(((j-1)*n_points_cyl)+i) cyl_peri_x(i+1) cyl_peri_y(i) cyl_peri_y(i+1) 1 2];
		
	end
	
	geom(n_points_cyl+var_len,:) = [2 cyl_cell{1,j}(n_points_cyl) cyl_cell{1,j}(1) cyl_cell{2,j}(n_points_cyl) cyl_cell{2,j}(1) 1 2];
	
	var_len = var_len + n_points_cyl;
	
end
%}

%% Creating CSG

% Rectangle
rect = [3 , 4 , -area_period/2 , area_period/2 , area_period/2 , -area_period/2 , -tot_height/2 , -tot_height/2 , tot_height/2 , tot_height/2];
ns = char('rect');
sf = 'rect';

geom = zeros(length(rect),n_cyl);

for i = 1:n_cyl

	geom(:,i) = [1 , cyl_cent_x(i) , cyl_cent_y(i) , r_cyl , zeros(1,length(rect) - 4)]';
	
	ns = char(ns,['circ',num2str(i)]);
	
	sf = [sf,'+',ns(i+1,:)];
	
end

ns = ns';

geom = [rect',geom];

%% Create Model, Geometry & Mesh
	
%{
% plot(points(1,:),points(2,:),'.')
% axis equal

% delaunayTriangulation

% geometryFromEdges(model,geom')

% pdegplot(model,'EdgeLabels','on')

% mesh = generateMesh(model,'Hmax',hmax,'Hgrad',1.05,'GeometricOrder','linear')

% [p,e,t] = initmesh(geom','hmax',hmax,'Hgrad',1.05,'MesherVersion','R2013a');

% h = pdemesh(p,e,t);
%}
[dl,bt] = decsg(geom,sf,ns);

% pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
% axis equal

model = createpde(1);

geometryFromEdges(model,dl);

mesh = generateMesh(model,'Hmax',hmax,'Hgrad',1.05);

% figure(2);pdeplot(model)

%% Triangle Manipulation

[p,e,t] = meshToPet(mesh);

tri_x_val = p(1,:);
tri_y_val = p(2,:);

n_nodes = length(tri_x_val);

n_tri = length(mesh.Elements(1,:));

B = [2 1 1 ; 1 2 1 ; 1 1 2]/12;

M = sparse(n_nodes,n_nodes);

% M = zeros(length(p(1,:)));

for i = 1:n_tri
	
	zone = t(end,i);
	
	if zone == 1
		
		diel_const = di_const1;
		
	else
		
		diel_const = di_const2;
		
	end
	
	xy_val = p(:,t(1:3,i));
	
	x1 = xy_val(1,1);
	x2 = xy_val(1,2);
	x3 = xy_val(1,3);
	y1 = xy_val(2,1);
	y2 = xy_val(2,2);
	y3 = xy_val(2,3);
	
	area_tri_k = abs((((x2 - x1) * (y3 - y1)) - ((y2 - y1) * (x3 - x1))))/2;
	
	du_dx = (y3 - y1)/((y3 - y1) * (x2 - y1) - (y2 - y1) * (x3 - y1));
	
	dv_dx = (y2 - y1)/((y2 - y1) * (x3 - y1) - (y3 - y1) * (x2 - y1));
	
	du_dy = -(x3 - x1)/((y3 - y1) * (x2 - y1) - (y2 - y1) * (x3 - y1));
	
	dv_dy = -(x2 - x1)/((y2 - y1) * (x3 - y1) - (y3 - y1) * (x2 - y1));
	
	d11=(-1*du_dx-1*dv_dx)*(-1*du_dx-1*dv_dx)+(-1*du_dy-1*dv_dy)*(-1*du_dy-1*dv_dy);
	d12=(-1*du_dx-1*dv_dx)*(+1*du_dx+0*dv_dx)+(-1*du_dy-1*dv_dy)*(+1*du_dy+0*dv_dy);
	d13=(-1*du_dx-1*dv_dx)*(+0*du_dx+1*dv_dx)+(-1*du_dy-1*dv_dy)*(+0*du_dy+1*dv_dy);   
	d21=(+1*du_dx+0*dv_dx)*(-1*du_dx-1*dv_dx)+(+1*du_dy+0*dv_dy)*(-1*du_dy-1*dv_dy);
	d22=(+1*du_dx+0*dv_dx)*(+1*du_dx+0*dv_dx)+(+1*du_dy+0*dv_dy)*(+1*du_dy+0*dv_dy);
	d23=(+1*du_dx+0*dv_dx)*(+0*du_dx+1*dv_dx)+(+1*du_dy+0*dv_dy)*(+0*du_dy+1*dv_dy);   
	d31=(+0*du_dx+1*dv_dx)*(-1*du_dx-1*dv_dx)+(+0*du_dy+1*dv_dy)*(-1*du_dy-1*dv_dy);
	d32=(+0*du_dx+1*dv_dx)*(+1*du_dx+0*dv_dx)+(+0*du_dy+1*dv_dy)*(+1*du_dy+0*dv_dy);
	d33=(+0*du_dx+1*dv_dx)*(+0*du_dx+1*dv_dx)+(+0*du_dy+1*dv_dy)*(+0*du_dy+1*dv_dy);   
	
	A = area_tri_k * [d11 d12 d13 ; d21 d22 d23 ; d31 d32 d33];
	
	Mk = (k0^2 * B - A/diel_const)*area_tri_k;
	
	M(i:i+2,i:i+2) = M(i:i+2,i:i+2) + Mk;
	
	i/n_tri*100
	
end








