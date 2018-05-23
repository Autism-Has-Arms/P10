classdef zone_determination < handle
	%
	
	properties
		
		PML
		env
		upper
		lower
		cyl
		cyl_whole
		cyl_split
		
	end
	
	
	methods
		
		function zone_det(obj,dl,bt,varargin)
			
			p = inputParser;
			
			default_values = {0 , 'circle'};
			valid_values = {[0,1] , {'circle','rectangle'}};
			check_values = {@(x) any(x == valid_values{1}) , @(x) any(strcmpi(x,valid_values{2}))};
			
			addParameter(p,'enable_surface',default_values{1},check_values{1});
			addParameter(p,'PML',default_values{1},check_values{1});
			addParameter(p,'ScatType',default_values{2},check_values{2});
			
			parse(p,varargin{:})
			
			%% Environment faces
			
			obj.env = find(sum(bt(:,1:(end - p.Results.enable_surface)),2) - p.Results.PML == 1).';
			
			if length(obj.env) ~= (1 + p.Results.enable_surface)
				
				warning('Could not determine faces of environment. Please insert them manually.')
				
				temp_fig = figure('Unit','normalized','Position',[(1-0.4-0.1) 0.25 0.4 0.5]);
				pdegplot(dl,'FaceLabels','on')
				axis equal
				
				envi = inputdlg({'Upper environment','Lower environment'},...
								 'Manual input',[1 30 ; 1 30]);
				
				obj.env = str2num(char(envi{:})).';
				
				close(temp_fig)
				
				key = [];
				
			elseif bt(obj.env(1),end)
				
				obj.env = fliplr(obj.env);
				
			end
			
			%% Cylinder faces
			
			cyl_geom = bt(:,(2 + p.Results.PML):(end - p.Results.enable_surface));
			
			if exist('key','var')
				
				cyl_geom(obj.env,:) = 0;
				
			end
			
			[row,col] = ind2sub(size(cyl_geom),find(cyl_geom));
			ind_cell = mat2cell(col == unique(col).',length(col),ones(1,length(unique(col))));
			
			for i = 1:length(ind_cell)
				
				cyl = row(ind_cell{i});
				
				if length(cyl) == 1
					
					obj.cyl_whole = [obj.cyl_whole cyl];
					
				else
					
					if bt(cyl(1),end)
						
						cyl = flipud(cyl);
						
					end
					
					obj.cyl_split = [obj.cyl_split cyl];
					
				end
				
			end
			
			obj.cyl = [obj.cyl_whole reshape(obj.cyl_split.',1,[])];
			
			%% Perfectly Matched Layer
			
			if p.Results.PML
			
				obj.PML = find(sum(bt(:,1:end - p.Results.enable_surface),2) == 1).';
				
				if bt(obj.PML(1),end)
					
					obj.PML = fliplr(obj.PML);
					
				end
				
			end
			
			if p.Results.enable_surface
				
				obj.upper = obj.env(bt(obj.env,end) == 0);
				obj.lower = obj.env(bt(obj.env,end) == 1);
				
				if p.Results.PML
					
					obj.upper = [obj.upper , obj.PML(bt(obj.PML,end) == 0)];
					obj.lower = [obj.lower , obj.PML(bt(obj.PML,end) == 1)];
					
				end
				
			end
			
		end
		
	end
	
end

