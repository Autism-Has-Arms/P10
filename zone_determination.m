classdef zone_determination < handle
	%
	
	properties
		
		PML
		env
		upper
		lower
		cyl
		
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
			
			obj.env = find(sum(bt(:,1:end - p.Results.enable_surface),2) - p.Results.PML == 1).';
			obj.cyl = find(sum(bt(:,1:end - p.Results.enable_surface),2) - p.Results.PML == 2).';
			
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

