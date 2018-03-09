clear all
close all
clc

%% Initialisation

if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2
	
	pct = disppct;
		
end


%% Parameters

di_const1 = 1;
di_const2 = 12;

hmax = 2 * pi * 10 / 20;

lambda = 700;

k0 = 2*pi/lambda;


%% Cylinder specifics

n_cyl = 7;
r_cyl = 10;
cyl_period = 50;


%% Area specifics

ul_spacing = 1500;
area_period = 30;
area_height = n_cyl * 300;
tot_height = ul_spacing + area_height;


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


%% Creating CSG

% Rectangle

rect = [3 , 4 , -area_period/2 , area_period/2 , area_period/2 , -area_period/2 , -tot_height/2 , -tot_height/2 , tot_height/2 , tot_height/2];
ns = char('rect');
sf = 'rect';

% Cylinders

geom = zeros(length(rect),n_cyl);

for i = 1:n_cyl

	geom(:,i) = [1 , cyl_cent_x(i) , cyl_cent_y(i) , r_cyl , zeros(1,length(rect) - 4)]';
	
	ns = char(ns,['circ',num2str(i)]);
	
	sf = [sf,'+',ns(i+1,:)];
	
end

ns = ns';

geom = [rect',geom];

%% Create Model, Geometry & Mesh

[dl,~] = decsg(geom,sf,ns);

% pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
% axis equal

model = createpde(1);

geometryFromEdges(model,dl);

mesh = generateMesh(model,'Hmax',hmax,'Hgrad',1.05,'GeometricOrder','linear');

% pdeplot(model);asd

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

bv = zeros(length(p),1);

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
	
	du_dx =  (y3 - y1)/((y3 - y1) * (x2 - y1) - (y2 - y1) * (x3 - y1));
	dv_dx =  (y2 - y1)/((y2 - y1) * (x3 - y1) - (y3 - y1) * (x2 - y1));
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
	
	Mk = k0^2 * B * area_tri_k - A/diel_const;
	
	
	% Triangles with a side touching top or bottom edge.
	
	if any(i == ind_top_edge) || any(i == ind_bot_edge)
		
		% x and y values of the three points in the i'th triangle.
		
		point_vec_ud = p(:,t(1:3,i));
		
		% Find which two indices have the same y value.
		
		ind_same_yval = logical(sum(repmat(point_vec_ud(2,:),3,1) == repmat(point_vec_ud(2,:)',1,3)) - 1);
		
		% Length is calculated as the difference between the corresponding
		% x values (since they have the same y value).
		
		edge_length = abs(diff(point_vec_ud(1,ind_same_yval)));
		
		% The y value is found by using the first index.
		
		y_val = point_vec_ud(2,find(ind_same_yval,1));
		
		H0 = exp(-1i * k0 * diel_const * y_val);
		
		bk = 1i * k0 * diel_const * H0 * edge_length;
		
		bv(t(ind_same_yval,i)) = bk;
		
		temp_mat = zeros(3);
		temp_mat(ind_same_yval,ind_same_yval) = [2 1 ; 1 2];
		
		C = 1i * edge_length * sqrt(diel_const) * k0 * (1/diel_const) * temp_mat/6;
		
		Mk = Mk + C;
		
	end
	
	M(t(1:3,i),t(1:3,i)) = M(t(1:3,i),t(1:3,i)) + Mk;
	
	
	
	if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2
	
		pct = disppct(i,n_tri,pct,1);
		
	else
		
		i/n_tri*100
	
	end
	
end



% Periodic boundary conditions

for i = 1:n_tri
	
	if any(i == ind_left_edge) || any(i == ind_right_edge)
		
		if any(i == ind_right_edge)
				
			ind_opposite = ind_left_edge;
				
		elseif any(i == ind_left_edge)
				
			ind_opposite = ind_right_edge;
				
		else
			
			error('i is neither in left or right.')
		
		end
		
		% x and y values of the three points in the i'th triangle.
		
		point_vec_lr = p(:,t(1:3,i));
		
		% Check which two indices in point_vec_lr have the same x value.
		
		ind_same_xval = logical(sum(repmat(point_vec_lr(1,:),3,1) == repmat(point_vec_lr(1,:)',1,3)) - 1);
		
		% The corresponding indices in the 'p' array.
		
		ind_in_p = t(ind_same_xval,i);
		
		% Compare with saved indices to see check if point has already been
		% considered.
		
		if ~isempty(ind_saved) && any(any((ind_in_p == ind_saved)'))
			
			% Removes an index if it has already been counted.
			
			ind_in_p = ind_in_p(not(any((ind_in_p == ind_saved)')));
			
		end
			
		% Check whether the exclusion of duplicate indices empties the
		% variable.
		
		if ~isempty(ind_in_p)
			
			continue
			
		end

		% Save the indices to not count them multiple times.

		ind_saved(length(ind_saved)+1:length(ind_saved)+length(ind_in_p)) = ind_in_p;

		% Corresponding y values on current edge.

		val_y_cur = p(2,ind_in_p);

		% y values of all points on the opposite side along their
		% indices in p. [y_value index_in_p].

		val_y_op = [p(2,t(1:3,ind_opposite)) ; reshape(t(1:3,ind_opposite),[1 numel(t(1:3,ind_opposite))])]';

		% Comparing y values and finding the indices of the minimum
		% values.

		[~,ind_min_op] = min(abs(val_y_cur - val_y_op(:,1)));

		% The indices in 'p' where these minimum values are located.
		% (The points on the opposite edge which are closest in height
		% to the selected ones on the current edge).

		ind_p_close = val_y_op(ind_min_op,2);

		% Change the corresponding row by introducing a row of zeros.

		M(ind_in_p,:) = 0;

		M(sub2ind(size(M),ind_in_p,ind_in_p)) = 1;

		M(sub2ind(size(M),ind_in_p,ind_p_close)) = -1;

		bv(ind_in_p) = 0;
		
	end
	
	if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2
	
		pct = disppct(i,n_tri,pct,2);
		
	else
		
		i/n_tri*100
	
	end
		
end

Hv = M\bv;

pdeplot(p,e,t,'xydata',abs(Hv))
axis equal
colormap gray


%% Indices of Triangles on Edges

function edge_index = edge_ind(mesh,x_or_y,num)
	
	% Gives indices of triangles on selected edge.

	switch x_or_y
		
		case 'x'
			
			x_or_y = 1;
			
		case 'y'
			
			x_or_y = 2;
			
		otherwise
			
			error("2nd input must be 'x' or 'y'.")
			
	end
	
	[p,e,t] = meshToPet(mesh);
	
	% Get the index in edge array 'e' where the x/y value first appears.
	
	ind_val = find(ismember([p(x_or_y,e(1,:));p(x_or_y,e(2,:))]',[num;num]','rows')',1);
	
	% Get the corresponding edge ID for the given x/y value.
	
	edge_id = e(5,ind_val);
	
	% All the edges with the given edge ID.
	
	apt_edges = (e(5,:) == edge_id);
	
	% Get the indices in 'p' of the 2 points which make up the edges.
	
	ind_edge_points = [e(1,apt_edges);e(2,apt_edges)];
	
	% Get indices of triangle array 't' which include two points on
	% selected edge.
	
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
	
end
