clear all
close all
clc


%% Initialisation.

r_i = load('Gold_refractive_index_file_J_C.m');
% r_i = load('Silver_refractive_index_file_J_C.m');

pct = disppct;

i_for = 0;

%% Parameters.

material_type = 'bulk';			% 'scatterers' or 'bulk'

scat_shape = 'circle';			% 'circle' or 'rectangle'.
placement_style = 'array';		% 'manual', 'array' or 'random'.

polarisation = 's';				% 'p' or 's'.

lambda = 700;

theta = 6*pi/4;

% Determines maximum size of elements. Therefore larger values of hmax
% creates fewer elements.
hmax = 4 .* 10;
hmin = hmax ./ 4;


r_env = 1500;					% Radius of environment.
r_scat = 12;					% Radius of scatterers.
r_PML = 800;					% Additional radius of PML. Total radius is r_env + r_PML.


% Surface parameters. If no surface is wanted, simply use the same
% refractive index as the environment.

n_surf = 1;						% Surface refractive index.
di_const_surf = n_surf.^2;
mag_const_surf = 1;
surface_y_coordinate = 0;		% Height of surface.

% Parameters for rectangular scatterers.
scat_width = 24;
scat_height = 24;
rounded_corners = true;

% Environment parameters.
n_env = 1;
di_const_env = n_env.^2;
mag_const_env = 1;

switch material_type

	case 'bulk'

		% Name of file containing effective parameters.			
		meta_str = ['s-pol,r_cyl=',num2str(r_scat),',wavelength=400-900,pattern=line.mat'];

		if exist(fullfile(pwd,meta_str),'file') ~= 2

			error(sprintf(['\nNo\n\n"',meta_str,'"\n\nin working directory.']))

		end

		meta_m = load(meta_str);

		di_const_scat = interp1(meta_m.wavelength,meta_m.permittivity,lambda);
		mag_const_scat = interp1(meta_m.wavelength,meta_m.permeability,lambda);

	case 'scatterers'

		n_scat = interp1(r_i(:,1),r_i(:,2)+r_i(:,3)*1i,lambda); % Cylinder
		di_const_scat = n_scat^2;
		mag_const_scat = 1;

end

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

theta_1 = theta;

if theta <= 3*pi/2

	theta_i = 3*pi/2 - theta_1;
	theta_t = asin(n_env/n_surf * sin(theta_i));
	theta_2 = 3*pi/2 - theta_t;

else

	theta_i = theta_1 - 3*pi/2;
	theta_t = asin(n_env/n_surf * sin(theta_i));
	theta_2 = theta_t + 3*pi/2;

end

refl = (n_surf * cos(theta_i) - n_env * cos(theta_t))/(n_env * cos(theta_t) + n_surf * cos(theta_i));
tran = (2 * n_surf * cos(theta_i))/(n_env * cos(theta_t) + n_surf * cos(theta_i));

k0 = 2*pi/lambda;

%% Centres of cylinders.

if strcmpi(placement_style,'array')

	obj_cent = cent_gen('MainStruct','circle','ScatStruct',scat_shape);

	obj_cent.gen_cent(rows_scat,period_scat,n_col,'StructShape','hexagonal','Scat_pm','+');

end


%% Creating CSG.

obj_csg = csg;

% Main structure

if strcmpi(material_type,'bulk')

	switch scat_shape

		case 'circle'

			scatterers = {[min(obj_cent.cent_x)-r_scat , max(obj_cent.cent_x)+r_scat] , [min(obj_cent.cent_y)-r_scat , max(obj_cent.cent_y)+r_scat]};

		case 'rectangle'

			scatterers = {[min(obj_cent.cent_x)-scat_width/2 max(obj_cent.cent_x)+scat_width/2] , [min(obj_cent.cent_y)-scat_height/2 max(obj_cent.cent_y)+scat_height/2]};

	end

	scat_width = 2.*scatterers{1}(2);

	scat_height = 2.*scatterers{2}(2);

	obj_cent.cent_x = 0;
	obj_cent.cent_y = 0;

	scat_shape = 'rectangle';

	placement_style = 'manual';

end

obj_csg.create_csg('circle',[0 0],r_env + r_PML);

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

obj_csg.create_csg('rectangle',2*(r_env + r_PML),[surface_y_coordinate , -(r_env + r_PML)]);


%% Create Model, Geometry & Mesh.

[dl,bt] = decsg(obj_csg.geom,obj_csg.sf,obj_csg.ns);

if strcmpi(scat_shape,'rectangle') && rounded_corners

	ind_edge_rem = zeros(1,4.*5.*size(obj_cent.cent_x,2));

	if strcmpi(material_type,'bulk')

		temp_dl = dl(8:9,:);

	else

		temp_dl = round(dl(8:9,:));

	end

	for j = 1:length(saved_corner_cent)

		corner_circ_cent = saved_corner_cent{j};
		circ_edge_keep = zeros(1,size(corner_circ_cent,2));

		for i = 1:size(corner_circ_cent,2)

			ind_circ = find(sum(temp_dl == corner_circ_cent(:,i)) == 2);

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

% figure
% pdegplot(dl,'FaceLabels','on','EdgeLabels','on')
% axis equal

model = createpde(1);

model.geometryFromEdges(dl);

mesh = generateMesh(model,'Hmax',hmax,'Hmin',hmin,'GeometricOrder','linear');


%% Zone determination.

obj_zone = zone_determination;

obj_zone.zone_det(dl,bt);

cyl_diel_const = di_const_scat * ones(1,length(obj_zone.scat));

env_diel_const = [di_const_env di_const_surf]; % Corresponding dielectric constant

cyl_mag_const = mag_const_scat * ones(1,length(obj_zone.scat));

env_mag_const = [mag_const_env mag_const_surf];

%% Triangle Manipulation.

[p,e,t] = meshToPet(mesh);

% figure
% pdemesh(p,e,t)

% Amount of times to refine the mesh in the scatterers.

n_refine = 2;

for i = 1:n_refine

	[p,e,t] = refinemesh(dl,p,e,t,obj_zone.scat);

end

% figure
% pdemesh(p,e,t)

n_tri = size(t,2);

B = [2 1 1 ; 1 2 1 ; 1 1 2]/12;

ind_PML = edge_ind({p,e,t},'r',r_env + r_PML);

n_for = 4;


%% Pre-calculation initialisations.

M = spalloc(size(p,2),size(p,2),3*length(p));

bv = zeros(length(p),1);


%% Calculations.

disp(['Triangle amount = ',num2str(length(t))])

i_for = i_for + 1;

for i = 1:n_tri

	zone = t(end,i);

	if any(zone == obj_zone.upper)

		diel_const = di_const_env;
		mag_const = mag_const_env;

	elseif any(zone == obj_zone.lower)

		diel_const = di_const_surf;
		mag_const = mag_const_surf;

	elseif any(zone == obj_zone.scat)

		diel_const = di_const_scat;
		mag_const = mag_const_scat;

	end

	ref_ind = sqrt(diel_const);

	xy_val = p(:,t(1:3,i));

	x1 = xy_val(1,1);
	x2 = xy_val(1,2);
	x3 = xy_val(1,3);
	y1 = xy_val(2,1);
	y2 = xy_val(2,2);
	y3 = xy_val(2,3);

	area_tri_i = abs((((x2 - x1) * (y3 - y1)) - ((y2 - y1) * (x3 - x1))))/2;

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

	fun_ang = cos(theta) * x_mid_tri + sin(theta) * y_mid_tri;

	E0 = exp(1i * k0 * ref_ind * fun_ang);

	r_0 = r_env;

	sigma_0 = 6*log(10)/(2*pi/lambda*r_PML.^3);

	if any(zone == obj_zone.PML)

		r_mid_tri = sqrt(x_mid_tri.^2 + y_mid_tri.^2);

		sigma = sigma_0/diel_const * (r_mid_tri - r_0).^2;

	else

		sigma = 0;

	end

	A = (area_tri_i * [d11 d12 d13 ; d21 d22 d23 ; d31 d32 d33]) ./ (1 + 1i * sigma).^2;

	switch polarisation

		case 'p'

			Mk = k0^2 * B * area_tri_i * mag_const - A/diel_const;

		case 's'

			Mk = k0^2 * B * area_tri_i * diel_const - A/mag_const;

	end

	M(t(1:3,i),t(1:3,i)) = M(t(1:3,i),t(1:3,i)) + Mk;

	% Progress of current for-loop.
	pct = disppct(i,n_tri,pct,i_for,n_for);

end

%% Perfectly Matched Layer.

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

	pct = disppct(i,n_tri,pct,i_for,n_for);

end


%% Contribution from all cylinder edges.

vec_normal = [];

i_for = i_for + 1;

for i = 1:length(e(1,:))

	if any(ismember(e(6:7,i),obj_zone.scat))

		l_or_r = logical(sum(e(6:7,i) == obj_zone.scat,2));

		vec = p(:,e(1:2,i));
		vec_orth = [vec(2,2) - vec(2,1) ; (vec(1,2) - vec(1,1))] .* (l_or_r - ~l_or_r);
		vec_normal = vec_orth/norm(vec_orth);

		if all(l_or_r)

			H_prime = (vec_normal(1)*cos(theta_2) + vec_normal(2)*sin(theta_2))...
					 * 1i * k0 * n_surf * tran * exp(1i*k0*n_surf*(cos(theta_2)*sum(vec(1,:))/2 + sin(theta_2)*sum(vec(2,:))/2));

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
					 * 1i * k0 * n_env * exp(1i*k0*n_env*(cos(theta_1)*sum(vec(1,:))/2 + sin(theta_1)*sum(vec(2,:))/2))...
					 + (vec_normal(1)*cos(theta_1) - vec_normal(2)*sin(theta_1))...
					 * 1i * k0 * n_env * refl * exp(1i*k0*n_env*(cos(theta_1)*sum(vec(1,:))/2 - sin(theta_1)*sum(vec(2,:))/2))...
					 ;...
					   (vec_normal(1)*cos(theta_2) + vec_normal(2)*sin(theta_2))...
					 * 1i * k0 * n_surf * tran * exp(1i*k0*n_surf*(cos(theta_2)*sum(vec(1,:))/2 + sin(theta_2)*sum(vec(2,:))/2))];

			switch polarisation

				case 'p'

					di_1 = env_diel_const(u_or_d);
					di_ref_1 = env_diel_const(1);
					di_2 = cyl_diel_const(l_or_r);
					di_ref_2 = env_diel_const(2);

					bk = (norm(vec_orth)/2) * (((1/di_1 - 1/di_ref_1) - (1/di_2 - 1/di_ref_2)) * H_prime(u_or_d));

				case 's'

					mag_1 = env_mag_const(u_or_d);
					mag_ref_1 = env_mag_const(1);
					mag_2 = cyl_mag_const(l_or_r);
					mag_ref_2 = env_mag_const(2);

					bk = (norm(vec_orth)/2) * ((1/mag_1 - 1/mag_ref_1) - (1/mag_2 - 1/mag_ref_2)) * H_prime(u_or_d);

			end

		end

		bv(e(1:2,i)) = bv(e(1:2,i)) + bk;

	end

	pct = disppct(i,length(e(1,:)),pct,i_for,n_for);

end


%% Contribution from triangles in all scatterers on the total field.

for i = 1:length(obj_zone.scat)

	tri_in_cyl = t(1:end-1,ismember(t(4,:),obj_zone.scat(i)));

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

			bk = tri_area * (2/6) .* (di_ref .* (1/cyl_diel_const(i) - 1/di_ref) - (cyl_mag_const(i) - mag_ref)) .* k0.^2 .* di_ref .* H0r;

		case 's'

			bk = tri_area * (2/6) .* (mag_ref .* (1/cyl_mag_const(i) - 1/mag_ref) - (cyl_diel_const(i) - di_ref)) .* k0.^2 .* mag_ref .* H0r;

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

		ref_ind = sqrt(di_const_env);

		H0v(i) = exp(1i * k0 * ref_ind * (cos(theta_1) * x_vec + sin(theta_1) * y_vec))...
			   + refl * exp(1i * k0 * ref_ind * (cos(theta_1) * x_vec - sin(theta_1) * y_vec));

	else

		ref_ind = sqrt(di_const_surf);

		H0v(i) = tran * exp(1i * k0 * ref_ind * (cos(theta_2) * x_vec + sin(theta_2) * y_vec));

	end

	pct = disppct(i,size(p,2),pct,i_for,n_for);

end

Hv = M\bv;

if false %n_env == n_surf	% Enable if surface line should be removed.

	temp_fig = figure('Unit','normalized','Position',[(1-0.4-0.1) 0.25 0.4 0.5]);
	pdegplot(dl,'EdgeLabels','on')
	axis equal

	edge_del = input('Surface edges to delete. Separate with comma:\n','s');

	close(temp_fig)

	[dl,bt] = csgdel(dl,bt,str2num(edge_del));

end

if ~exist('var_object','var')

	figure('Unit','normalized','Position',[0.05 0.25 0.4 0.5])
	pdeplot(p,e,t,'xydata',abs(Hv.'))%,'Zdata',abs(Hv))
	colormap parula
	hold on
	pdegplot(dl)
	axis equal

	Hvn = Hv + H0v;

	figure('Unit','normalized','Position',[(1-0.4-0.05) 0.25 0.4 0.5])
	pdeplot(p,e,t,'xydata',abs(Hvn))%,'Zdata',abs(Hv))
	colormap parula
	hold on
	pdegplot(dl)
	axis equal

end

%% Plotting values of line down through structure.

angle_peri = linspace(0.01,2*pi,10000);

line_x = cos(angle_peri) * r_env;

line_y = sin(angle_peri) * r_env;

int_F = pdeInterpolant(p,t,Hv);

line_abs = (abs(evaluate(int_F,[line_x;line_y])).^2)*(r_env - 100);

switch polarisation

	case 'p'

		line_abs = line_abs .* n_env;
		line_abs(line_y < surface_y_coordinate) = line_abs(line_y < surface_y_coordinate) .* n_surf ./ n_env;

	case 's'

		line_abs = line_abs ./ n_env;
		line_abs(line_y < surface_y_coordinate) = line_abs(line_y < surface_y_coordinate) .* n_env ./ n_surf;

end

curve_area = trapz(angle_peri,line_abs/r_env);
	
figure
plot(angle_peri,line_abs)
% scatter3(cos(angle_peri)*r_env,sin(angle_peri)*r_env,line_abs,1,line_abs)


