clear all
close all
clc


%% Initialisation

if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2

	pct = disppct;

end

var_object = 'hmax = linspace(3,15,200)';

if exist('var_object','var')
	
	var_string = [extractBefore(var_object,' = ') , ' = var_array(k);'];
	
	var_array = eval(extractAfter(var_object,' = '));
	
	var_len = length(var_array);
	
else
	
	var_len = 1;
	
end

reflectance = zeros(1,length(var_len));

transmittance = zeros(1,length(var_len));

distance = zeros(1,length(var_len));

ref_index = zeros(1,length(var_len));

wavelength = zeros(1,length(var_len));

cyl_amount = zeros(1,length(var_len));

cyl_radius = zeros(1,length(var_len));

cyl_periods = zeros(1,length(var_len));

for k = 1:var_len

	%% Parameters

	lambda = 700;
	
	r_i = load('Gold_refractive_index_file_J_C.m');
	% r_i = load('Silver_refractive_index_file_J_C.m');
	
	% Determines maximum size of elements. Therefore larger values of hmax
	% creates fewer elements. 
	hmax = 2 * pi * 10 / 25;
	
	% Cylinder specifics

	n_cyl = 10;
	r_cyl = 10;
	cyl_period = 30;
	
	% Area specifics

	ul_spacing = 1400;
	area_width = 30;
	
% 	hmax = var_array(k);
	
	if exist('var_string','var')
		
		eval(var_string)
		
	end

	% Minor calculations
	
	n1 = 1;
	n2 = interp1(r_i(:,1),r_i(:,2)+r_i(:,3)*1i,lambda);
	di_const1 = n1^2;
	di_const2 = n2^2;

	% Test of script duartion:
	% Halving hmax from 2pi*10/10 to 2pi*10/20
	% increased the scrpt duration 4 times, while halving it again from
	% 2pi*10/20 to 2pi*10/40 increased the duration 10 times.

	k0 = 2*pi/lambda;


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

	cent_cyl = [cyl_cent_x ; cyl_cent_y];

	area_height = 2*(max(cyl_cent_y) + cyl_period/2) + 100;
	tot_height = ul_spacing + area_height;


	%% Creating CSG

	a = csg;

	% Rectangle

	a.create_csg('rectangle',area_width,tot_height);

	% Gold rectangle

	% gold_rect = [3 , 4 , area_width/2 , -area_width/2 , -area_width/2 , area_width/2 , area_height/2 , area_height/2 , -area_height/2 , -area_height/2];
	% ns = char(ns,'rect_g');

	% Cylinders

	for i = 1:n_cyl

		a.create_csg('circle',cent_cyl(:,i),r_cyl);

	end
	

	%% Create Model, Geometry & Mesh

	[dl,~] = decsg(a.geom,a.sf,a.ns);

	% pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
	% axis equal

	model = createpde(1);

	geometryFromEdges(model,dl);

	mesh = generateMesh(model,'Hmax',hmax,'Hgrad',1.05,'GeometricOrder','linear');

	% pdeplot(model);

	%% Triangle Manipulation

	[p,e,t] = meshToPet(mesh);

	n_tri = length(mesh.Elements(1,:));

	B = [2 1 1 ; 1 2 1 ; 1 1 2]/12;

	M = sparse(size(mesh.Nodes,2),size(mesh.Nodes,2));

	ind_top_edge = edge_ind(mesh,'y',tot_height/2);

	ind_bot_edge = edge_ind(mesh,'y',-tot_height/2);

	ind_right_edge = edge_ind(mesh,'x',area_width/2);

	ind_left_edge = edge_ind(mesh,'x',-area_width/2);

	ind_saved = [];

	bv = zeros(length(p),1);

	%% Calculations

	for i = 1:n_tri

		zone = t(end,i);
		
		if zone == 2 %|| zone == 1

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

		du_dx =  (y3 - y1)/((y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1));
		dv_dx =  (y2 - y1)/((y2 - y1) * (x3 - x1) - (y3 - y1) * (x2 - x1));
		du_dy = -(x3 - x1)/((y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1));
		dv_dy = -(x2 - x1)/((y2 - y1) * (x3 - x1) - (y3 - y1) * (x2 - x1));

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

			if any(i == ind_top_edge)

				% The y value is found by using the first index.

				y_val = point_vec_ud(2,find(ind_same_yval,1));

				H0 = exp(-1i * k0 * sqrt(diel_const) * y_val);

				bk = 1i * k0 * diel_const * H0 * edge_length;

				bv(t(ind_same_yval,i)) = bv(t(ind_same_yval,i)) + bk; %<-- Husk vinkelafhængig.

			end

			temp_mat = zeros(3);
			temp_mat(ind_same_yval,ind_same_yval) = [2 1 ; 1 2];

			C = 1i * edge_length * sqrt(diel_const) * k0 * (1/diel_const) * temp_mat/6;

			Mk = Mk + C;

		end

		M(t(1:3,i),t(1:3,i)) = M(t(1:3,i),t(1:3,i)) + Mk;

		if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2

			pct = disppct(i,n_tri,pct,2*k-1,2*var_len);

		else

			i/n_tri*100

		end

	end


	% Periodic boundary conditions

	for i = 1:n_tri

		if any(i == ind_left_edge)

			ind_opposite = ind_right_edge;

			% x and y values of the three points in the i'th triangle.

			point_vec_lr = p(:,t(1:3,i));

			% Check which two indices in point_vec_lr have the same x value.

			ind_same_xval = logical(sum(repmat(point_vec_lr(1,:),3,1) == repmat(point_vec_lr(1,:)',1,3)) - 1);

			% The corresponding indices in the 'p' array.

			ind_in_p = t(ind_same_xval,i);

			% Compare with saved indices to see check if point has already been
			% considered.

			if ~isempty(ind_saved) && any(any((ind_in_p == repmat(ind_saved,length(ind_in_p),1))'))

				% Removes an index if it has already been counted.

				ind_in_p = ind_in_p(not(any((ind_in_p == repmat(ind_saved,length(ind_in_p),1))')));

			end

			% Check whether the exclusion of duplicate indices empties the
			% variable. If true, skip iteration.

			if isempty(ind_in_p)

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


		end

		if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2

			pct = disppct(i,n_tri,pct,2*k,2*var_len);

		else

			i/n_tri*100

		end

	end

	% Hv = lsqminnorm(M,bv);

	Hv = M\bv;

	% pdeplot(p,e,t,'xydata',abs(Hv))
	% hold on
	% pdegplot(dl)
	% axis equal
	% colormap gray 


	%% Plotting values of line down through structure

	line_x = linspace(0,0,500);

	line_y = linspace(-tot_height,tot_height,500);

	% Interpolant

	int_F = pdeInterpolant(p,t,Hv);

	line_abs = abs(evaluate(int_F,[line_x;line_y]));

	% figure(2)

	% plot(line_y,line_abs)

	% axis([-tot_height tot_height -1.5 1.5])


	%% Calculating transmittance and reflectance etc.

	val_y_t = evaluate(int_F,[0 ; -tot_height/2]);

	val_y_r = evaluate(int_F,[0 ; tot_height/2]);

	y0 = area_height/2;

	H0 = exp(-1i*k0*n1*tot_height/2);

	reflectance(k) = (val_y_r - H0) * exp(2i*k0*n1*y0) * exp(-1i*k0*n1*tot_height/2);

	transmittance(k) = val_y_t * exp(1i*k0*n1*y0) * exp(1i*k0*n1*(tot_height/2 + y0));

	distance(k) = 2 * y0;

	ref_index(k) = n2;

	wavelength(k) = lambda;
	
	cyl_amount(k) = n_cyl;
	
	cyl_radius(k) = r_cyl;
	
	cyl_periods(k) = cyl_period;
	
end

%% Writing to file

save(['Var=' , extractBefore(var_object,' = ') , ',[' , num2str(min(var_array)) , ...
	',' , num2str(max(var_array)) , '],n=' , num2str(var_len) , '.mat'],...
	'wavelength','ref_index','distance','reflectance','transmittance','var_array',...
	'cyl_periods','area_width','cyl_amount','cyl_radius');

%{

header = {'Wavelength' 'Refractive index' 'Distance' 'Reflectance' 'Transmittance'};
comma_header = [header;repmat({','},1,numel(header))];
comma_header = comma_header(:)';
text_header = cell2mat(comma_header);
csvv = [lambda,n2,2*(max(cyl_cent_y) + cyl_period/2),refl,val_y_t];
%}

dispstat('Finished.','keepprev','timestamp');

%% Indices of Triangles on Edges

function edge_index = edge_ind(mesh,x_or_y,num)
	
	% Gives indices of triangles (in array 't') on selected edge.

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
	
	ind_val = find(ismember(round([p(x_or_y,e(1,:));p(x_or_y,e(2,:))]',5),[num;num]','rows')',1);
	
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