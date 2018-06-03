clear all
close all
clc


%% Initialisation

r_i = load('Gold_refractive_index_file_J_C.m');
% r_i = load('Silver_refractive_index_file_J_C.m');

if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2

	pct = disppct;

end

var_object = 'lambda = linspace(300,400,2)';

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

scat_cross = zeros(1,length(var_len));

i_for = 0;

for k = 1:var_len

	%% Parameters
	
	main_structure = 'circle';	% 'circle' or 'rectangle'.
	scat_shape = 'circle';			% 'circle' or 'rectangle'.
	material_type = 'scatterers';			% 'scatterers' or 'bulk'
	
	polarisation = 's';				% 'p' or 's'.
	
	PML = true;						% Include Perfectly-Matched Layer.
	
	enable_surface = true;			% Enable surface.
	
	geometric_order = 'linear';		% 'linear' or 'quadratic'.

	lambda = 700;
	
	theta = 6*pi/4;
	
	% Determines maximum size of elements. Therefore larger values of hmax
	% creates fewer elements.
	hmax = 3;
	hmin = hmax ./ 4;
	
	n1 = 1;
	n2 = interp1(r_i(:,1),r_i(:,2)+r_i(:,3)*1i,lambda); % Cylinder
	
	if strcmpi(scat_shape,'circle')
	
		r_scat = 12; %sqrt((24.^2 - 4.^2 .* (4 - pi))./pi);	% Radius of scatterers.
		
	elseif strcmpi(scat_shape,'rectangle')

		scat_width = 24;
		scat_height = 24;
		rounded_corners = true;

	end
	
	switch material_type
		
		case 'bulk'
			
			meta_str = ['s-pol,r_cyl=',num2str(r_scat),',wavelength=400-900,pattern=line.mat'];
			
			if exist(fullfile(pwd,meta_str),'file') ~= 2
				
				error(sprintf(['\nNo\n\n"',meta_str,'"\n\nin working directory.']))
				
			end
			meta_m = load(meta_str);
			di_const1 = interp1(meta_m.wavelength,meta_m.permittivity,lambda);
			di_const3 = n2^2;
			mag_const1 = interp1(meta_m.wavelength,meta_m.permeability,lambda);
			mag_const3 = 1;
			
		case 'scatterers'
			
			di_const1 = n1^2;
			di_const3 = n2^2;
			mag_const1 = 1;
			mag_const3 = 1;
			
	end
	
	if strcmp(main_structure,'circle')
		
		hmax = hmax * 10;
		hmin = hmin * 10;
		r_env = 1500;					% Radius of environment.
		
		placement_style = 'manual';		% 'manual', 'array' or 'random'.
		
		if strcmpi(placement_style,'manual')
			
			scatterers = {12,0,0};
			
		elseif strcmpi(placement_style,'array')
			
			n_col = 13;
			rows_scat = 13;
			period_scat = 30;
			
		elseif strcmpi(placement_style,'random')
			
			cyl_n = 100;
			cyl_dist = 2.*r_scat + 4;
			max_tries = 100000;
			
			
			area_x = [-100 , 100];
			area_y = [-100 , 100];
			
			cent_x = zeros(1,cyl_n);
			cent_y = zeros(1,cyl_n);
			
			i = 1;
			while i <= cyl_n
				
				cent_x(i) = area_x(1) + (area_x(2) - area_x(1)).*rand;
				cent_y(i) = sin(acos(cent_x(i)./area_x(2))) .* (area_y(1) + (area_y(2) - area_y(1)).*rand);
				
				inter_dist = [cent_x(1:(i-1)) - cent_x(i) ; cent_y(1:(i-1)) - cent_y(i)];
				
				if i == 1 || all(sqrt(inter_dist(1,:).^2 + inter_dist(2,:).^2) > cyl_dist)
					
					i = i + 1;
					
					count = 0;
					
				else
					
					count = count + 1;
					
					if count == max_tries
						
						cent_x(cent_x == 0) = [];
						cent_x(end) = [];
						cent_y(cent_y == 0) = [];
						cent_y(end) = [];
						
						break
						
					end
					
				end
				
			end
			
			obj_cent = cent_gen;
			obj_cent.cent_x = cent_x;
			obj_cent.cent_y = cent_y;
			
		else
			
			if strcmpi(scat_shape,'rectangle')

				scatterers = {scat_width scat_height};	% {Width , Height} or {[x_start x_end] , [y_start , y_end]} or another permutation.

			elseif strcmpi(scat_shape,'circle')

				scatterers = {r_scat 0 0};	% {Radii , Centre x-coordinates , Centre y_coordinates] of scatterer.

			end
			
		end
		
		if PML
			
			r_PML = 800;	% Distanace from r_circ.
			
		end
		
		if enable_surface
			
			surface_y_coordinate = 0;
			
			n2 = 1.5;	% Glass.
			di_const2 = n2^2;
			mag_const2 = 1;
			
			theta_1 = theta;
			
			if theta <= 3*pi/2
				
				theta_i = 3*pi/2 - theta_1;
				theta_t = asin(n1/n2 * sin(theta_i));
				theta_2 = 3*pi/2 - theta_t;
				
			else
				
				theta_i = theta_1 - 3*pi/2;
				theta_t = asin(n1/n2 * sin(theta_i));
				theta_2 = theta_t + 3*pi/2;
				
			end
			
			refl = (n2 * cos(theta_i) - n1 * cos(theta_t))/(n1 * cos(theta_t) + n2 * cos(theta_i));
			tran = (2 * n2 * cos(theta_i))/(n1 * cos(theta_t) + n2 * cos(theta_i));
			
		end
		
	elseif strcmp(main_structure,'rectangle')
		
		placement_style = 'array';
		
		rows_scat = 11;
		period_scat = 30;
		
		ul_spacing = 1400;
		area_width = 4*period_scat;
		
		cyl_pattern = 'hexagonal'; % ['line','hexagonal'].

		if strcmp(cyl_pattern,'hexagonal') && r_scat >= (period_scat/sqrt(2))

			error(['Overlapping cylinders. Radius (r_cyl) must be below r = ' , num2str(round(period_scat/sqrt(2),2)) , '. (Hexagonal structure)'])

		elseif strcmp(cyl_pattern,'line') && r_scat >= (period_scat/2)

			error(['Overlapping cylinders. Radius (r_cyl) must be below r = ' , num2str(round(period_scat/2,2)) , '. (Line structure)'])
			
		end
		
	end
	
	% Minor calculations
	
	if exist('var_string','var')
		
		eval(var_string)
		
	end

	% Test of script duartion:
	% Halving hmax from 2pi*10/10 to 2pi*10/20
	% increased the scrpt duration 4 times, while halving it again from
	% 2pi*10/20 to 2pi*10/40 increased the duration 10 times.

	k0 = 2*pi/lambda;


	if length(lambda) == 1 || k == 1	
		%% Centres of cylinders

		if strcmpi(placement_style,'array')
			
			obj_cent = cent_gen('MainStruct',main_structure,'ScatStruct',scat_shape);
		
			if strcmpi(obj_cent.main_struct,'rectangle')

				if rows_scat == 0

					tot_height = ul_spacing;

				else

					obj_cent.gen_cent(rows_scat,period_scat,area_width,'StructShape',cyl_pattern);

					area_height = 2*max(obj_cent.cent_y) + period_scat;

					tot_height = ul_spacing + area_height;

				end

			else

				obj_cent = cent_gen;

				obj_cent.gen_cent(rows_scat,period_scat,n_col,'StructShape','hexagonal','Scat_pm','-');
			
			end
			
		end


		%% Creating CSG

		obj_csg = csg;

		% Main structure
		
		if strcmp(main_structure,'rectangle')
		
			obj_csg.create_csg('rectangle',area_width,tot_height);
			obj_csg.sf = 'rect';
			
			% Cylinders

			for i = 1:length(obj_cent.cent_y)

	% 			a.create_csg('rectangle',area_width,[r_cyl , cent_cyl(2,i)]);

				obj_csg.create_csg('circle',[obj_cent.cent_x(i) ; obj_cent.cent_y(i)],r_scat);

			end
		
		elseif strcmp(main_structure,'circle')
			
			if PML
		
				obj_csg.create_csg('circle',[0 0],r_env + r_PML);
				
			end
			
			obj_csg.sf = 'circ';
			
			obj_csg.create_csg('circle',[0 , 0],r_env);
			
			if strcmpi(placement_style,'array') || strcmpi(placement_style,'random')

				for i = 1:length(obj_cent.cent_x)
					
					if strcmpi(scat_shape,'circle')

						obj_csg.create_csg('circle',[obj_cent.cent_x(i) , obj_cent.cent_y(i)],r_scat);
						
					else
						
						obj_csg.create_csg('rectangle',scat_width,scat_height,'Centre',[obj_cent.cent_x(i) , obj_cent.cent_y(i)]);
						
					end

				end
				
			else
			
				for i = 1:size(scatterers,1)
					
					if strcmpi(scat_shape,'rectangle')
						
						obj_csg.create_csg('rectangle',scatterers{i,1},scatterers{i,2});
						
					elseif strcmpi(scat_shape,'circle')

						obj_csg.create_csg('circle',[scatterers{i,2} , scatterers{i,3}],scatterers{i,1});
						
					end

				end
				
			end
			
			if strcmpi(scat_shape,'rectangle') && rounded_corners

				r_corner_circ = (scat_width+scat_height)/12;
				
				x_offset = scat_width/2 - r_corner_circ;
				y_offset = scat_height/2 - r_corner_circ;
				
				saved_corner_cent = cell(length(obj_cent.cent_x),1);

				for i = 1:length(obj_cent.cent_x)

					corner_circ_cent = [ obj_cent.cent_x(i) + x_offset .* [1 -1 -1 1] ; obj_cent.cent_y(i) + y_offset .* [1 1 -1 -1] ];

					obj_csg.create_csg('circle',[corner_circ_cent(1,1) , corner_circ_cent(2,1)],r_corner_circ)
					obj_csg.create_csg('circle',[corner_circ_cent(1,2) , corner_circ_cent(2,2)],r_corner_circ)
					obj_csg.create_csg('circle',[corner_circ_cent(1,3) , corner_circ_cent(2,3)],r_corner_circ)
					obj_csg.create_csg('circle',[corner_circ_cent(1,4) , corner_circ_cent(2,4)],r_corner_circ)
					
					saved_corner_cent{i} = corner_circ_cent;

				end

			end
			
% 			obj_csg.create_csg('circle',[0 , 0],165);

% 			obj_csg.create_csg('circle',[15 , 15],r_scat);

% 			obj_csg.create_csg('circle',[0 , 0],r_scat * (1 + 1/4));

% 			obj_csg.create_csg('circle',[0 , 0],r_circ * (1 - 1/30));

			if enable_surface
				
				if PML
				
					obj_csg.create_csg('rectangle',2*(r_env + r_PML),[surface_y_coordinate , -(r_env + r_PML)]);
					
				else
					
					obj_csg.create_csg('rectangle',2*r_env,[surface_y_coordinate , -r_env]);
					
				end
				
			end
			
		end


		%% Create Model, Geometry & Mesh

		[dl,bt] = decsg(obj_csg.geom,obj_csg.sf,obj_csg.ns);

		if strcmpi(scat_shape,'rectangle') && rounded_corners

			ind_edge_rem = zeros(1,4.*5.*size(obj_cent.cent_x,2));
			
			for j = 1:length(saved_corner_cent)
				
				corner_circ_cent = saved_corner_cent{j};
				circ_edge_keep = zeros(1,size(corner_circ_cent,2));
			
				for i = 1:size(corner_circ_cent,2)

					ind_circ = find(sum(round(dl(8:9,:)) == corner_circ_cent(:,i)) == 2);

					face_lr = dl(6:7,ind_circ);

					face_not_in_circ = ~[all(all(repmat(face_lr(1,:),size(corner_circ_cent,2),1) == repmat(face_lr(1,:),size(corner_circ_cent,2),1).')) ; ...
										 all(all(repmat(face_lr(2,:),size(corner_circ_cent,2),1) == repmat(face_lr(2,:),size(corner_circ_cent,2),1).'))];

					face_pos = sum(repmat(face_lr(face_not_in_circ,:),size(corner_circ_cent,2),1) == repmat(face_lr(face_not_in_circ,:),size(corner_circ_cent,2),1).') == 1;

					face_label = face_lr(face_not_in_circ,face_pos);

					ind_of_face = logical(sum(face_lr == face_label));

					circ_edge_keep = ind_circ(ind_of_face);

					ind_edge_face = find(sum(dl(6:7,:) == face_label));

					ind_edge_rem((j-1)*5*4 + 5*(i-1) + (1:2)) = ind_edge_face(ind_edge_face ~= circ_edge_keep);

					ind_edge_rem((j-1)*5*4 + 5*(i-1) + (3:5)) = ind_circ(~ind_of_face);

				end
				
			end
			
			[dl,bt] = csgdel(dl,bt,ind_edge_rem);
		
		end

% 		figure
% 		pdegplot(dl)%,'FaceLabels','on','EdgeLabels','on')
% 		axis equal

		model = createpde(1);

		model.geometryFromEdges(dl);

		mesh = generateMesh(model,'Hmax',hmax,'Hmin',hmin,'GeometricOrder',geometric_order);
		

		%% Zone determination
		
		switch main_structure
			
			case 'circle'
				
				obj_zone = zone_determination;

				obj_zone.zone_det(dl,bt,'enable_surface',enable_surface,'PML',PML);
			
			case 'rectangle'
				
				obj_zone = zone_determination;
				
				unique_zones = unique(dl(7,:)); 
   
				unique_zones = unique_zones(~unique_zones == 0); 

				[~ , zone_majority_ind] = max(histc(dl(7,:),unique_zones)); 

				obj_zone.env = unique_zones(zone_majority_ind);
				
				zone_scat = 1:size(bt,2);
				obj_zone.cyl = zone_scat(zone_scat ~= obj_zone.env);
  
		end
		
		cyl_diel_const = di_const3 * ones(1,length(obj_zone.cyl));
		
		env_diel_const = di_const1; % Corresponding dielectric constant
		
		cyl_mag_const = mag_const3 * ones(1,length(obj_zone.cyl));
		
		env_mag_const = mag_const1;
		
		if enable_surface
			
			env_diel_const = [env_diel_const di_const2];
			
			env_mag_const = [env_mag_const mag_const2];
			
		end
		
		
		%% Triangle Manipulation

		[p,e,t] = meshToPet(mesh);
		
% 		figure
% 		pdemesh(p,e,t)
		
% 		[p,e,t] = MidPointFix({p,e,t});
		
		n_refine = 2;

		for i = 1:n_refine

			[p,e,t] = refinemesh(dl,p,e,t,obj_zone.cyl);
		
		end
		
% 		figure
% 		pdemesh(p,e,t)
		
		n_tri = size(t,2);

		B = [2 1 1 ; 1 2 1 ; 1 1 2]/12;
		
		if strcmp(main_structure,'rectangle')

			ind_top_edge = edge_ind({p,e,t},'y',tot_height/2);

			ind_bot_edge = edge_ind({p,e,t},'y',-tot_height/2);

			ind_right_edge = edge_ind({p,e,t},'x',area_width/2);

			ind_left_edge = edge_ind({p,e,t},'x',-area_width/2);
			
			n_for = 2*var_len;
			
		elseif strcmp(main_structure,'circle')
			
			if PML
				
				ind_PML = edge_ind({p,e,t},'r',r_env + r_PML);
				
				n_for = 2*var_len;
				
			else
			
				ind_peri_edge = edge_ind({p,e,t},'r',r_env);

				n_for = var_len;
			
			end
			
			if enable_surface
				
				n_for = n_for + 2*var_len;
				
			end
			
		end
		
	end
	
	
	%% Pre-calculation initialisations
	
	M = spalloc(size(p,2),size(p,2),3*length(p));

	bv = zeros(length(p),1);

	
	%% Calculations
	
	if var_len == 1
		
		disp(['Triangle amount = ',num2str(length(t))])
		
	end

	i_for = i_for + 1;
	
	for i = 1:n_tri

		zone = t(end,i);
		
		if enable_surface
			
			if any(zone == obj_zone.upper)
				
				diel_const = di_const1;
				mag_const = mag_const1;
			
			elseif any(zone == obj_zone.lower)
				
				diel_const = di_const2;
				mag_const = mag_const2;
				
			elseif any(zone == obj_zone.cyl)
				
				diel_const = di_const3;
				mag_const = mag_const3;
				
			end
			
		else
		
			if any(zone == obj_zone.env)

				diel_const = di_const1;
				mag_const = mag_const1;

			else

				diel_const = di_const3;
				mag_const = mag_const3;

			end
			
		end
		
		ref_ind = sqrt(diel_const);
		
		if strcmp(mesh.GeometricOrder,'linear')

			xy_val = p(:,t(1:3,i));
			
		else
			
			xy_val = p(:,t(1:6,i));
			
		end

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
		
		% x and y coordinates of middle of triangle.
		x_mid_tri = ((max(xy_val(1,:)) - min(xy_val(1,:)))/2) + min(xy_val(1,:));
		y_mid_tri = ((max(xy_val(2,:)) - min(xy_val(2,:)))/2) + min(xy_val(2,:));
		
		r = sqrt(x_mid_tri.^2 + y_mid_tri.^2);
		
		fun_ang = cos(theta) * x_mid_tri + sin(theta) * y_mid_tri;

		E0 = exp(1i * k0 * ref_ind * fun_ang);
			
		A = area_tri_k * [d11 d12 d13 ; d21 d22 d23 ; d31 d32 d33];
		
		if PML
			
			r_0 = r_env;
			
			sigma_0 = 6*log(10)/(2*pi/lambda*r_PML^3);%4*log(10)/(k0*(r_PML.^2)*(diel_const.^2));

			if any(t(end,i) == obj_zone.PML)
			
				sigma = sigma_0/diel_const * (r - r_0).^2;
				
			else
				
				sigma = 0;
				
			end

			A = A / (1 + 1i * sigma).^2;
			
		else
			
			switch polarisation
				
				case 'p'
			
					bk = - k0.^2 * (1/diel_const - 1/di_const1) * E0 * 2 * area_tri_k / 6;
					
				case 's'
					
					bk = - k0.^2 * (diel_const - di_const1) * E0 * 2 * area_tri_k / 6;
					
			end
			
			bv(t(1:3,i)) = bv(t(1:3,i)) + bk.';

		end

		switch polarisation
			
			case 'p'
			
				Mk = k0^2 * B * area_tri_k * mag_const - A/diel_const;
			
			case 's'
			
				Mk = k0^2 * B * area_tri_k * diel_const - A/mag_const;
			
		end


		%% Triangles with a side touching top or bottom edge

		if strcmp(main_structure,'rectangle') && (any(i == ind_top_edge) || any(i == ind_bot_edge))

			% Find which two indices have the same y value.

			ind_same_yval = logical(sum(repmat(xy_val(2,:),3,1) == repmat(xy_val(2,:),3,1)') - 1);

			% Length is calculated as the difference between the corresponding
			% x values (since they have the same y value).

			edge_length = abs(diff(xy_val(1,ind_same_yval)));

			% The x and y values are found.

			x_val = xy_val(1,ind_same_yval);
			y_val = xy_val(2,ind_same_yval);

			k_x = k0 * ref_ind * cos(theta);
			k_y = k0 * ref_ind * sin(theta);
			
			H0 = exp(-1i * k_y * y_val) .* exp(1i * k_x * x_val);
			
% 			if (pi < theta && theta <= 2*pi) && any(i == ind_top_edge)

				bk = 1i * k_y * ref_ind * H0 * edge_length;
				
% 			end
			
% 			if  (0 <= theta && theta < pi) && any(i == ind_bot_edge)

% 				bk = 1i * k_y * ref_ind * H0 * edge_length;
				
% 			end
			
			bv(t(ind_same_yval,i)) = bv(t(ind_same_yval,i)) + bk.'; %<-- Husk vinkelafhængig.

			temp_mat = zeros(3);
			temp_mat(ind_same_yval,ind_same_yval) = [2 1 ; 1 2];

			C = 1i * edge_length * sqrt(diel_const) * k0 * (1/diel_const) * temp_mat/6;

			Mk = Mk + C;

		end
		
		
		%% Circle geometry
		
		if strcmp(main_structure,'circle') && ~PML && any(i == ind_peri_edge)
			
			% Find which two indices are on the periphery.
			
			% Overall length is found.
			ind_peri = sqrt(xy_val(1,:).^2 + xy_val(2,:).^2);
			
			% Corner-corner average length.
			cc_length = sqrt(diff(xy_val(1,[1 2 3 1])).^2 + diff(xy_val(2,[1 2 3 1])).^2);
			cc_length = sum(cc_length)/length(cc_length);
			
			% Round result to 2 decimal places.
			ind_peri = abs(round(ind_peri,2) - r_env) < cc_length/20;
			
			% Distance between points is found.
			inter_point_dist = diff(xy_val(:,ind_peri),1,2);

			% Length is calculated by using the Pythagorean Theorem.
			edge_length = min(sqrt(sum(inter_point_dist.^2)));

			% The y value is found by using the first index.

			fun_ang = cos(theta) * xy_val(1,ind_peri) + sin(theta) * xy_val(2,ind_peri);

			E0 = exp(1i * k0 * sqrt(diel_const) * fun_ang);
				
			bk = ((1 - fun_ang/r_env) * 1i * k0 * sqrt(diel_const) - 1/(2 * r_env)) .* E0 * edge_length/(2 * diel_const);

			bk = [2 * bk(1) + bk(2) , bk(1) + 2 * bk(2)]/3;

			bv(t(ind_peri,i)) = bv(t(ind_peri,i)) + bk.';

			temp_mat = zeros(3);
			temp_mat(ind_peri,ind_peri) = [2 1 ; 1 2];

			C = (edge_length/6 * (1i * k0 * sqrt(diel_const) - 1/(2 * r_env)) * temp_mat)/diel_const;

			Mk = Mk + C;

		end

		M(t(1:3,i),t(1:3,i)) = M(t(1:3,i),t(1:3,i)) + Mk;

		if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2

			pct = disppct(i,n_tri,pct,i_for,n_for);

		else

			i/n_tri*100

		end

	end


	%% Periodic boundary conditions

	if strcmp(main_structure,'rectangle')
		
		ind_saved = [];
		
		i_for = i_for + 1;
		
		for i = 1:n_tri

			if  any(i == ind_left_edge)
				
				x_val = p(1,t(1:3,i));

				ind_opposite = ind_right_edge;

				% Check which two indices in xy_val have the same x value.
				ind_same_xval = logical(sum(repmat(x_val,3,1) == repmat(x_val',1,3)) - 1);

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
				ind_p_close = val_y_op(ind_min_op,2);
				
				k_x = k0 * diel_const.^2 * sin(theta);
 
				% Setting a row to zero.
				M(ind_in_p,:) = 0;

				M(sub2ind(size(M),ind_in_p,ind_in_p)) = 1;
				
				M(sub2ind(size(M),ind_in_p,ind_p_close)) = -exp(1i*k_x*area_width);

				bv(ind_in_p) = 0;


			end

			if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2

				pct = disppct(i,n_tri,pct,i_for,n_for);

			else

				i/n_tri*100

			end

		end
		
	end
	
	%% Perfectly Matched Layer
	
	if strcmp(main_structure,'circle') && PML
		
		ind_saved = [];
		
		i_for = i_for + 1;
		
		for i = 1:n_tri
			
			if any(i == ind_PML)
				
				xy_val = p(:,t(1:3,i));
				
				% Find which two indices are on the periphery.
			
				% Overall length is found.
				ind_peri = sqrt(xy_val(1,:).^2 + xy_val(2,:).^2);

				% Corner-corner average length.
				cc_length = sqrt(diff(xy_val(1,[1 2 3 1])).^2 + diff(xy_val(2,[1 2 3 1])).^2);
				cc_length = sum(cc_length)/length(cc_length);

				% Find which two of the indices are on the periphery.
				ind_peri = abs(round(ind_peri,2) - (r_env + r_PML)) < cc_length/20;
				
				% Indices in 'p'.
				ind_in_p = t(ind_peri,i);
				
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
				
				bv(ind_in_p) = 0;
				
				M(ind_in_p,:) = 0;
				M(sub2ind(size(M),ind_in_p,ind_in_p)) = 1;
				
			end
			
			if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2

				pct = disppct(i,n_tri,pct,i_for,n_for);

			else

				i/n_tri*100

			end
			
		end
		
	end
	

	%% Contribution from all cylinder edges.
	
	if enable_surface
	
		vec_normal = [];
		
% 		j = 0;
		i_for = i_for + 1;

		for i = 1:length(e(1,:))

			if any(ismember(e(6:7,i),obj_zone.cyl))

% 				j = j + 1;
				l_or_r = logical(sum(e(6:7,i) == obj_zone.cyl,2));

				vec = p(:,e(1:2,i));
				vec_orth = [vec(2,2) - vec(2,1) ; (vec(1,2) - vec(1,1))] .* (l_or_r - ~l_or_r);
				vec_normal = vec_orth/norm(vec_orth);
% 				normal_vec(:,j) = vec_normal+mean(vec,2);
				
				if all(l_or_r)
					
					H_prime = (vec_normal(1)*cos(theta_2) + vec_normal(2)*sin(theta_2))...
							 * 1i * k0 * n2 * tran * exp(1i*k0*n2*(cos(theta_2)*sum(vec(1,:))/2 + sin(theta_2)*sum(vec(2,:))/2));
					
					if vec_normal(2) == 1
						
						u_or_d = logical([1 ; 0]);
						
					else
						
						u_or_d = logical([0 ; 1]);
						
					end
					
					switch polarisation

						case 'p'
							
							di_1 = cyl_diel_const(u_or_d);
							di_ref_1 = env_diel_const(u_or_d);
							di_2 = cyl_diel_const(~u_or_d);
							di_ref_2 = env_diel_const(~u_or_d);

							bk = (norm(vec_orth)/2)*((1/di_1 - 1/di_ref_1)-(1/di_2 - 1/di_ref_2))*H_prime;

						case 's'
							
							mag_1 = cyl_mag_const(u_or_d);
							mag_ref_1 = env_mag_const(u_or_d);
							mag_2 = cyl_mag_const(~u_or_d);
							mag_ref_2 = env_mag_const(~u_or_d);

							bk = (norm(vec_orth)/2)*((1/mag_1 - 1/mag_ref_1)-(1/mag_2 - 1/mag_ref_2))*H_prime;

					end
					
				else

					u_or_d = logical(sum(e(6:7,i) == obj_zone.env,2));

					H_prime = [(vec_normal(1)*cos(theta_1) + vec_normal(2)*sin(theta_1))...
							 * 1i * k0 * n1 * exp(1i*k0*n1*(cos(theta_1)*sum(vec(1,:))/2 + sin(theta_1)*sum(vec(2,:))/2))...
							 + (vec_normal(1)*cos(theta_1) - vec_normal(2)*sin(theta_1))...
							 * 1i * k0 * n1 * refl * exp(1i*k0*n1*(cos(theta_1)*sum(vec(1,:))/2 - sin(theta_1)*sum(vec(2,:))/2))...
							 ;...
							   (vec_normal(1)*cos(theta_2) + vec_normal(2)*sin(theta_2))...
							 * 1i * k0 * n2 * tran * exp(1i*k0*n2*(cos(theta_2)*sum(vec(1,:))/2 + sin(theta_2)*sum(vec(2,:))/2))];

					switch polarisation

						case 'p'
							
							di_1 = env_diel_const(u_or_d);
							di_ref_1 = env_diel_const(1);
							di_2 = cyl_diel_const(l_or_r);
							di_ref_2 = env_diel_const(2);

							bk = (norm(vec_orth)/2)*(((1/di_1 - 1/di_ref_1) - (1/di_2 - 1/di_ref_2))*H_prime(u_or_d));

						case 's'
							
							mag_1 = env_mag_const(u_or_d);
							mag_ref_1 = env_mag_const(1);
							mag_2 = cyl_mag_const(l_or_r);
							mag_ref_2 = env_mag_const(2);

							bk = (norm(vec_orth)/2)*((1/mag_1 - 1/mag_ref_1) - (1/mag_2 - 1/mag_ref_2))*H_prime(u_or_d);

					end
					
				end
				
				bv(e(1:2,i)) = bv(e(1:2,i)) + bk;

			end
			
			if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2

				pct = disppct(i,length(e(1,:)),pct,i_for,n_for);

			else

				i/n_tri*100

			end

		end
	
	%% Contribution from triangles in all scatterer on the total field.
	
		for i = 1:length(obj_zone.cyl)
			
			tri_in_cyl = t(1:end-1,ismember(t(4,:),obj_zone.cyl(i)));

			tri_cyl_x = p(1,tri_in_cyl);
			tri_cyl_y = p(2,tri_in_cyl);
			xy_123_vals = [1:3:length(tri_cyl_x) ; 2:3:length(tri_cyl_x) ; 3:3:length(tri_cyl_x)];
			tri_area = abs((tri_cyl_y(xy_123_vals(3,:)) - tri_cyl_y(xy_123_vals(1,:)))...
						.* (tri_cyl_x(xy_123_vals(2,:)) - tri_cyl_x(xy_123_vals(1,:)))...
						 - (tri_cyl_y(xy_123_vals(2,:)) - tri_cyl_y(xy_123_vals(1,:)))...
						.* (tri_cyl_x(xy_123_vals(3,:)) - tri_cyl_x(xy_123_vals(1,:))))/2;

			x_avg = (tri_cyl_x(xy_123_vals(1,:)) + tri_cyl_x(xy_123_vals(2,:)) + tri_cyl_x(xy_123_vals(3,:)))/3;
			y_avg = (tri_cyl_y(xy_123_vals(1,:)) + tri_cyl_y(xy_123_vals(2,:)) + tri_cyl_y(xy_123_vals(3,:)))/3;

			if mean(y_avg) > surface_y_coordinate

				di_ref = env_diel_const(1);

				mag_ref = env_mag_const(1);

				H0r = exp(1i*k0*(sqrt(di_ref))*(cos(theta_1)*x_avg + sin(theta_1)*y_avg))...
					+ refl*exp(1i*k0*(sqrt(di_ref))*(cos(theta_1)*x_avg - sin(theta_1)*y_avg));

			else

				di_ref = env_diel_const(2);

				mag_ref = env_mag_const(2);

				H0r = tran*exp(1i*k0*(sqrt(di_ref))*(cos(theta_2)*x_avg + sin(theta_2)*y_avg));

			end

			switch polarisation

				case 'p'

					bk = tri_area*(2/6) .* (di_ref .*(1/cyl_diel_const(i) - 1/di_ref) - (cyl_mag_const(i) - mag_ref)) .* k0.^2 .* di_ref .* H0r;

				case 's'

					bk = tri_area*(2/6) .* (mag_ref .*(1/cyl_mag_const(i) - 1/mag_ref) - (cyl_diel_const(i) - di_ref)) .* k0.^2 .* mag_ref .* H0r;

			end
			
			for j = 1:length(bk)

				bv(tri_in_cyl(:,j)) = bv(tri_in_cyl(:,j)) + bk(j);
				
			end

		end
		
		H0v = zeros(size(p,2),1);
		
		i_for = i_for + 1;
		
		for i = 1:size(p,2)
			
			x_vec = p(1,i);
			y_vec = p(2,i);
			
			if y_vec >= surface_y_coordinate
				
				ref_ind = sqrt(di_const1);
				
				H0v(i) = exp(1i * k0 * ref_ind * (cos(theta_1) * x_vec + sin(theta_1) * y_vec))...
					   + refl * exp(1i * k0 * ref_ind * (cos(theta_1) * x_vec - sin(theta_1) * y_vec));
				
			else
				
				ref_ind = sqrt(di_const2);
				
				H0v(i) = tran * exp(1i * k0 * ref_ind * (cos(theta_2) * x_vec + sin(theta_2) * y_vec));
				
			end
			
			if exist('disppct.m','file') == 2 && exist('dispstat.m','file') == 2

				pct = disppct(i,size(p,2),pct,i_for,n_for);

			else

				i/n_tri*100

			end
			
		end
	
	end
	
	Hv = M\bv;
	
	if var_len == 1
	
		figure('Unit','normalized','Position',[0.05 0.25 0.4 0.5])
		pdeplot(p,e,t,'xydata',abs(Hv.'))%,'Zdata',abs(Hv))
		colormap gray %parula
		hold on
		pdegplot(dl)
		axis equal
	
	if strcmpi(main_structure,'rectangle')
		
		H0vxy = cos(theta) * p(1,:) + sin(theta) * p(2,:);
		H0v = exp(1i*k0*n1*H0vxy).';
		
	end
		
	Hvn = Hv + H0v;

	figure('Unit','normalized','Position',[(1-0.4-0.05) 0.25 0.4 0.5])
	pdeplot(p,e,t,'xydata',abs(Hvn))%,'Zdata',abs(Hv))
	colormap gray %parula
	hold on
	pdegplot(dl)
	axis equal
	
	end

% 	caxis([0 1.5])

% 	figure;pdemesh(p,e,t)
% 	hold on
% 	plot(normal_vec(1,:),normal_vec(2,:),'*')

	%% Plotting values of line down through structure
	
	angle_peri = linspace(0,2*pi,10000);
	
	line_x = cos(angle_peri) * r_env;
	
	line_y = sin(angle_peri) * r_env;
	
	%{
	line_x = linspace(0,0,500);

	line_y = linspace(-tot_height,tot_height,500);
	%}
	% Interpolant
	
% 	fun_ang = cos(theta) .* line_x + sin(theta) .* line_y;

% 	E0_line = exp(1i * k0 * sqrt(1) * line_x).'; % fun_ang <=> line_x

	int_F = pdeInterpolant(p,t,Hv);

% 	line_abs = abs(evaluate(int_F,[line_x;line_y]) - E0_line).^2;
	
	line_abs = (abs(evaluate(int_F,[line_x;line_y])).^2)*(r_env - 100);
	
	switch polarisation
		
		case 'p'
	
			line_abs = line_abs .* n1;
			
		case 's'
			
			line_abs = line_abs ./ n1;
			
	end
	
	if enable_surface
		
		switch polarisation
			
			case 'p'
	
				line_abs(line_y < surface_y_coordinate) = line_abs(line_y < surface_y_coordinate) .* n2 ./ n1;
				
			case 's'
				
				line_abs(line_y < surface_y_coordinate) = line_abs(line_y < surface_y_coordinate) .* n1 ./ n2;
				
		end
		
	end
	
% 	line_abs = line_abs / max(line_abs);
	
	curve_area = trapz(angle_peri,line_abs/r_env);
	
% 	figure

% 	scatter3(cos(angle_peri)*r_circ,sin(angle_peri)*r_circ,line_abs,1,line_abs)
% 	plot(angle_peri,line_abs)
% 	asd
% 	plot(line_y,line_abs)

% 	axis([-tot_height tot_height -1.5 1.5])
	scat_cross(k) = curve_area;
	
end

plot(var_array,scat_cross)

%{
	%% Calculating transmittance and reflectance etc.

	val_y_t = evaluate(int_F,[0 ; -tot_height/2]);

	val_y_r = evaluate(int_F,[0 ; tot_height/2]);

	y0 = area_height/2;

	H0 = exp(-1i*k0*n1*tot_height/2);

	reflectance(k) = (val_y_r - H0) * exp(2i*k0*n1*y0) * exp(-1i*k0*n1*tot_height/2);

	transmittance(k) = val_y_t * exp(1i*k0*n1*y0) * exp(1i*k0*n1*(tot_height/2 - y0));

	distance(k) = 2 * y0;

	ref_index(k) = n2;

	wavelength(k) = lambda;
	
	cyl_amount(k) = rows_cyl;
	
	cyl_radius(k) = r_circ;
	
	cyl_periods(k) = cyl_period;
	
	scat_cross(k) = curve_area;
	
end
%}
%{
%% Writing to file

if logical(exist('var_array','var'))

	save(['Var=' , extractBefore(var_object,' = ') , ',[' , num2str(min(var_array)) , ...
	',' , num2str(max(var_array)) , '],n=' , num2str(var_len) , '.mat'],...
	'wavelength','ref_index','distance','reflectance','transmittance','var_array',...
	'cyl_periods','area_width','cyl_amount','cyl_radius');

else
	
	save(['Var=none,n=' , num2str(var_len) , '.mat'],...
	'wavelength','ref_index','distance','reflectance','transmittance',...
	'cyl_periods','area_width','cyl_amount','cyl_radius');

end

%{

header = {'Wavelength' 'Refractive index' 'Distance' 'Reflectance' 'Transmittance'};
comma_header = [header;repmat({','},1,numel(header))];
comma_header = comma_header(:)';
text_header = cell2mat(comma_header);
csvv = [lambda,n2,2*(max(cyl_cent_y) + cyl_period/2),refl,val_y_t];
%}

dispstat('Finished.','keepprev','timestamp');
%}
%% Indices of Triangles on Edges

function edge_index = edge_ind(pet,x_or_y,num)
	
	% Gives indices of triangles (in array 't') on selected edge.
	
	[p,e,t] = pet{:};
	
	switch x_or_y
		
		case 'x'
			
			point_val = [p(1,e(1,:)) ; p(1,e(2,:))];
			
		case 'y'
			
			point_val = [p(2,e(1,:)) ; p(2,e(2,:))];
			
		case 'r'
			
			point_val = sqrt([p(1,e(1,:)) ; p(1,e(2,:))].^2 + [p(2,e(1,:)) ; p(2,e(2,:))].^2);
			
		otherwise
			
			error("2nd input must be 'x', 'y' or 'r'.")
			
	end
	
	% Get the corresponding edge ID for the given x/y/r value.
	
	edge_id = unique(e(5,sum(abs(point_val - num) < 1) == 2));
	
	% All the edges with the given edge ID.
	
	apt_edges = logical(sum(e(5,:) == edge_id',1));
	
	% Get the indices in 'p' of the 2 points which make up the edges.
	
	ind_edge_points = [e(1,apt_edges);e(2,apt_edges)];
	
	% Choosing only the 't' indices where two points appear.
	
	edge_index = find(sum(ismember(t(1:3,:),ind_edge_points)) == 2);
	
	%{
	% Plot edge points
	xv = p(1,e(1:2,apt_edges));
	yv = p(2,e(1:2,apt_edges));
	hold on
	plot(xv,yv,'black*')
	%}
	
end

%% Fixing midpoint in quadratic

function [p,e,t] = MidPointFix(pet)

	[p,e,t] = pet{:};
	
	% The two edge points and the point in-between.
	
	edge_points = [1 2 4 ; 2 3 5 ; 3 1 6];
	
	for i = 1:size(edge_points,1)
	
		point_val_x = reshape(p(1,t(edge_points(i,:),:).'),[],3);
		point_val_y = reshape(p(2,t(edge_points(i,:),:).'),[],3);
		
		% Intermediate values
		
		int_val_x = (point_val_x(:,2) - point_val_x(:,1))/2 + point_val_x(:,1);
		int_val_y = (point_val_y(:,2) - point_val_y(:,1))/2 + point_val_y(:,1);
		
		p(1,t(edge_points(i,3),diff([point_val_x(:,3),int_val_x].') ~= 0)) = int_val_x(diff([point_val_x(:,3),int_val_x].') ~= 0);
		p(2,t(edge_points(i,3),diff([point_val_y(:,3),int_val_y].') ~= 0)) = int_val_y(diff([point_val_y(:,3),int_val_y].') ~= 0);
		
	end

end
