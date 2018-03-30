classdef cent_gen
	%CENT_GEN Generates the x and y coordinates of cylinder centres.
	
	properties
		
		cent_x
		cent_y
		
	end
	
	methods
		
		function obj = cent_gen(rows_cyl,cyl_period,area_width,struct_shape)
			
			insert = @(num, arr, pos) cat(2,  arr(1:pos-1), num, arr(pos:end));
			
			if mod(rows_cyl,2) == 0				% Checks if n_cyl is an even number.

				cyl_cent_y = zeros(1,rows_cyl);

				for i = 1:rows_cyl/2

					cyl_cent_y(2*i-1) = i * cyl_period - cyl_period/2;
					cyl_cent_y(2*i) = - i * cyl_period + cyl_period/2;

				end

			else
				
				cyl_cent_y = zeros(1,rows_cyl);

				for i = 1:(rows_cyl-1)/2
					
					cyl_cent_y(2*i) = i * cyl_period;
					cyl_cent_y((2*i)+1) = - i * cyl_period;
					
				end

			end
			
			if strcmp(struct_shape,'line')
			
				cyl_cent_x = zeros(1,rows_cyl);
				
			else
				
				switch mod(rows_cyl,4)
					
					case 0
						
						tot_cyl = rows_cyl + 2*floor(rows_cyl/4);
					
					case 1
						
						tot_cyl = rows_cyl + 2*floor(rows_cyl/4);
						
					case 2
						
						tot_cyl = rows_cyl + 2*floor(rows_cyl/4) + 1;
						
					case 3
						
						tot_cyl = rows_cyl + 2*floor(rows_cyl/4) + 2;
						
				end
			
				cyl_cent_x = zeros(1,tot_cyl);

				count = 1;

				for i = 1:rows_cyl

					if any(mod(i,4) == [0 , 1])

						cyl_cent_x(count) = 0;

						count = count+1;

					elseif any(mod(i,4) == [2 , 3])

						cyl_cent_x([count,count+1]) = [area_width/2 , -area_width/2];

						cyl_cent_y = insert(cyl_cent_y(count),cyl_cent_y,count+1);

						count = count+2;

					else

						error("Mod stuff doesn't work")

					end

				end
				
			end
			
			obj.cent_x = cyl_cent_x;
			obj.cent_y = cyl_cent_y;
			
		end
		
	end
	
end

