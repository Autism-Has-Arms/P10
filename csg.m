classdef csg < handle
	%CSG Creates Constructive Solid Geometry.
	%   Detailed explanation goes here
	
	properties
		
		geom
		ns
		sf = 'rect'
		
	end
	
	methods
		
		function create_csg(obj,type,a,b)
			
			%CSG Construct an instance of this class
			%   Detailed explanation goes here
			
			switch type
				
				case 'rectangle'
					
					type_char = 'rect';
					
					if length(b) == 1
					
						temp_geom = [3 , 4 , a/2 , -a/2 , -a/2 , a/2 , b/2 , b/2 , -b/2 , -b/2];
						
					else
						
						temp_geom = [3 , 4 , a/2 , -a/2 , -a/2 , a/2 , (b(1)/2)+b(2) , (b(1)/2)+b(2) , (-b(1)/2)+b(2) , (-b(1)/2)+b(2)];
						
					end
					
					obj.geom = [obj.geom' ; temp_geom]';
					
				case 'circle'
					
					type_char = 'circ';
					
					temp_geom = [1 , a(1) , a(2) , b , zeros(1,6)];
					obj.geom = [obj.geom' ; temp_geom]';
					
			end
			
			if isempty(obj.ns)
				
				obj.ns = char(type_char)';
				obj.sf = obj.ns';
				
			else
				
				temp_str = [type_char,num2str( length( obj.ns(1,:) ))];
				obj.ns = char(obj.ns',temp_str)';
				
			end
		
		end
		
	end
	
end

