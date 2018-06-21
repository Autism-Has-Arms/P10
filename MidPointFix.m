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