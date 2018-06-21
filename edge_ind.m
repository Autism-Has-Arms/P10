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