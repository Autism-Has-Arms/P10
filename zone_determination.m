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
		
		function zone_det(obj,bt,varargin)
			
			p = inputParser;
			
			default_values = 0;
			valid_values = [0,1];
			check_values = @(x) any(x == valid_values);
			
			addParameter(p,'enable_surface',default_values,check_values);
			addParameter(p,'PML',default_values,check_values);
			
			parse(p,varargin{:})
			
			%% Environment faces
			
			obj.env = find(sum(bt(:,1:(end - p.Results.enable_surface)),2) - p.Results.PML == 1).';
			
			%% Cylinder faces
			
			cyl_geom = bt(:,(2 + p.Results.PML):(end - p.Results.enable_surface));
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
			
			if p.Results.PML
			
				obj.PML = find(sum(bt(:,1:end - p.Results.enable_surface),2) == 1).';
				
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

