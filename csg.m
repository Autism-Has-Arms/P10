classdef csg < handle
	% csg Creates Constructive Solid Geometry (CSG).
	%	Class must be called to create object and gain access 
	%	to containing methods/functions.
	%
	% csg Properties:
	%	geom - The geometry of the structure in csg format.
	%	ns - The custom name of the structure.
	%	sf - The formula of the geometries.
	%
	% csg Methods:
	%	create_csg - Geometric Generation.
	%	Every call to this function appends data to the object's
	%	properties.
	%	If object type is a rectangle, (a) is width and (b) is
	%	height.
	%	If object type is a circle, (a) is centre coordinates and
	%	(b)	is radius.
	
	properties
		
		geom
		ns
		sf
		
	end
	
	methods
		
		function create_csg(obj,type,a,b)
			
			% create_csg - Geometric Generation.
			
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

