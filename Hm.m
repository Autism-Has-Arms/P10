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
% n_points_ud = 50;
% n_points_lr = floor((ul_spacing + area_height)/10);


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

% [p,e,t] = initmesh(geom','hmax',hmax,'Hgrad',1.05,'MesherVersion','R2013a');

% h = pdemesh(p,e,t);
%}

[dl,bt] = decsg(geom,sf,ns);

pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
% axis equal

model = createpde(1);

geometryFromEdges(model,dl);

mesh = generateMesh(model,'Hmax',hmax,'Hgrad',1.05,'GeometricOrder','linear');

pdeplot(model)

%% Triangle Manipulation

[p,e,t] = meshToPet(mesh);

point_x_val = p(1,:);
point_y_val = p(2,:);

n_nodes = length(point_x_val);

n_tri = length(mesh.Elements(1,:));

B = [2 1 1 ; 1 2 1 ; 1 1 2]/12;

M = sparse(n_nodes,n_nodes);

ind_top_edge = edge_ind(mesh,'y',tot_height/2);

ind_bot_edge = edge_ind(mesh,'y',-tot_height/2);

ind_right_edge = edge_ind(mesh,'x',area_period/2);

ind_left_edge = edge_ind(mesh,'x',-area_period/2);

ind_saved = [];

%% Calculations

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
	
	d11 = (-1 * du_dx - 1 * dv_dx) * (-1 * du_dx - 1 * dv_dx) + (-1 * du_dy - 1 * dv_dy) * (-1 * du_dy - 1 * dv_dy);
	d12 = (-1 * du_dx - 1 * dv_dx) * (+1 * du_dx + 0 * dv_dx) + (-1 * du_dy - 1 * dv_dy) * (+1 * du_dy + 0 * dv_dy);
	d13 = (-1 * du_dx - 1 * dv_dx) * (+0 * du_dx + 1 * dv_dx) + (-1 * du_dy - 1 * dv_dy) * (+0 * du_dy + 1 * dv_dy);   
	d21 = (+1 * du_dx + 0 * dv_dx) * (-1 * du_dx - 1 * dv_dx) + (+1 * du_dy + 0 * dv_dy) * (-1 * du_dy - 1 * dv_dy);
	d22 = (+1 * du_dx + 0 * dv_dx) * (+1 * du_dx + 0 * dv_dx) + (+1 * du_dy + 0 * dv_dy) * (+1 * du_dy + 0 * dv_dy);
	d23 = (+1 * du_dx + 0 * dv_dx) * (+0 * du_dx + 1 * dv_dx) + (+1 * du_dy + 0 * dv_dy) * (+0 * du_dy + 1 * dv_dy);
	d31 = (+0 * du_dx + 1 * dv_dx) * (-1 * du_dx - 1 * dv_dx) + (+0 * du_dy + 1 * dv_dy) * (-1 * du_dy - 1 * dv_dy);
	d32 = (+0 * du_dx + 1 * dv_dx) * (+1 * du_dx + 0 * dv_dx) + (+0 * du_dy + 1 * dv_dy) * (+1 * du_dy + 0 * dv_dy);
	d33 = (+0 * du_dx + 1 * dv_dx) * (+0 * du_dx + 1 * dv_dx) + (+0 * du_dy + 1 * dv_dy) * (+0 * du_dy + 1 * dv_dy);
	
	A = area_tri_k * [d11 d12 d13 ; d21 d22 d23 ; d31 d32 d33];
	
	Mk = (k0^2 * B - A/diel_const) * area_tri_k;
	
	% Top or bottom triangles.
	
	if any(i == ind_top_edge) || any(i == ind_bot_edge)
		
		point_vec = p(1:2,t(1:3,i));
		
		if point_vec(2,1) == point_vec(2,2)
			
			edge_length = abs(point_vec(1,1) - point_vec(1,2));
			
			y_val = point_vec(2,1);
			
		elseif point_vec(2,1) == point_vec(2,3)
			
			edge_length = abs(point_vec(1,1) - point_vec(1,3));
			
			y_val = point_vec(2,1);
			
		elseif point_vec(2,2) == point_vec(2,3)
			
			edge_length = abs(point_vec(1,2) - point_vec(1,3));
			
			y_val = point_vec(2,2);
			
		else
			
			error('No ifs.')
		
		end
		
		H0 = exp(-1i * k0 * diel_const * y_val);
		
		bk = 1i * k0 * diel_const * H0 * edge_length;
		
		C = 1i * edge_length * sqrt(diel_const) * k0 * (1/diel_const) * [2 1 0 ; 1 2 0 ; 0 0 0]/6;
		
		Mk = Mk + C;
		
	end
	
	% Periodic boundary conditions
	
	if any(i == ind_left_edge) || any(i == ind_right_edge)
		
		if any(i == ind_right_edge)
				
				ind_opposite = ind_left_edge;
				
		elseif any(i == ind_left_edge)
				
				ind_opposite = ind_right_edge;
				
		else
			
			error('i is neither in left or right.')
		
		end
		
		% x and y values of the three points in the i'th triangle.
		
		point_vec = p(1:2,t(1:3,i));
		
		% Check which two indices in point_vec have the same x value.
		
		ind_same_xval = logical(sum(point_vec(1,:) == point_vec(1,:)') - 1);
		
		% The corresponding indices in the 'p' array.
		
		ind_in_p = t(ind_same_xval,i);
		
		% Save these to not count them multiple times.
		
		ind_saved(length(ind_saved)+1:length(ind_saved)+length(ind_in_p)) = ind_in_p;
		
		% Corresponding y values.
		
		val_y_cor = p(2,ind_in_p);
		
		% Consider only opposite points with same x value.
		
		
		
		% Find points on opposite side with the same y value.
		
		%ind_op_side = 
		
	end
	
	M(t(1:3,i),t(1:3,i)) = M(t(1:3,i),t(1:3,i)) + Mk;
	
	i/n_tri*100
	
end


%% Edges

function edge_index = edge_ind(mesh,x_or_y,num)
	
	% Gives indices of points on selected edge.

	switch x_or_y
		
		case 'x'
			
			x_or_y = 1;
			
		case 'y'
			
			x_or_y = 2;
			
		otherwise
			
			error("3rd input must be 'x' or 'y'.")
			
	end
	
	[p,e,t] = meshToPet(mesh);
	
	% Get the index in edge array 'e' where the x/y value first appears.
	
	ind_val = find(ismember([p(x_or_y,e(1,:));p(x_or_y,e(2,:))]',[num;num]','rows')',1);
	
	% Get the corresponding edge ID for the given x/y value.
	
	edge_id = e(5,ind_val);
	
	% All the edges with the given edge ID.
	
	apt_edges = e(5,:) == edge_id;
	
	% Get the indices of the 2 points which make up the edges.
	
	ind_edge_points = [e(1,apt_edges);e(2,apt_edges)];
	
	% Get indices of triangle array 't' which include two points on
	% selected edge 
	
	ind_tri = ceil(find(ismember(t(1:3,:),ind_edge_points))/3);
	
	uni_ind = unique(ind_tri);
	
	ind_no = [uni_ind,histc(ind_tri,uni_ind)];
	
	% Choosing only the indices which appear twice.
	
	edge_index = uni_ind(ind_no(:,2) == 2)';
	
	%{
	% Plot edge points
	xv = p(1,e(1:2,apt_edges));
	yv = p(2,e(1:2,apt_edges));
	hold on
	plot(xv,yv,'r*')
	%}
	
	%{
	
	point_x_val = p(1,:);
	point_y_val = p(2,:);
	
	if edge == 'x'
		
		point_val = point_x_val;
		
	elseif edge == 'y'
		
		point_val = point_y_val;
		
	else
		
		error("x or y edges?")
		
	end

	% find((e(6,:) == 0) + (e(7,:) == 0) == 1)

	% Index in e where edge interfaces with subdomain 0 (exterior).
	e_bool_ex = (e(6,:) == 0) + (e(7,:) == 0) == 1;

	% Finding the corresponding indices in p for first and second point of edge.
	point_1_ind = e(1,e_bool_ex);
	point_2_ind = e(2,e_bool_ex);

	% Inserting indices in point_val to find values of points on edge.
	edge_point_1_val = point_val(point_1_ind);
	edge_point_2_val = point_val(point_2_ind);
	
	% Indices where point 1 and 2 have the same x/y value & has edge value.
	same_val_ind_1 = point_1_ind((abs(edge_point_1_val - edge_point_2_val) < 1) & (abs(edge_point_1_val - num) < 1));
	same_val_ind_2 = point_2_ind((abs(edge_point_1_val - edge_point_2_val) < 1) & (abs(edge_point_1_val - num) < 1));

% 	% The edge of box is chosen.
% 
% 	point_1_ext_ind = same_val_ind_1(point_val(same_val_ind_1) == num)
% 	point_2_ext_ind = same_val_ind_2(point_val(same_val_ind_2) == num)

	counter = 0;
	
	tri_up_bound = 0;
	
	for i = 1:length(t(1,:))
		
			if (ismember(t(1,i),same_val_ind_1) && ismember(t(2,i),same_val_ind_2)) || (ismember(t(1,i),same_val_ind_1) && ismember(t(3,i),same_val_ind_2)) || (ismember(t(2,i),same_val_ind_1) && ismember(t(3,i),same_val_ind_2))
				
				
				counter = counter + 1;

				tri_up_bound(counter) = i;

			end

	end
	
	edge_index = tri_up_bound;
	
	%}

end

