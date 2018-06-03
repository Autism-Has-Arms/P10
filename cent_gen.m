classdef cent_gen < handle
	%CENT_GEN Generates the x and y coordinates of cylinder centres.
	
	properties
		
		main_struct
		scat_struct
		cent_x
		cent_y
		
	end
	
	methods
		
		function obj = cent_gen(varargin)
			
			p = inputParser;
			
			default_values = 'circle';
			valid_values = {'circle','rectangle'};
			check_values = @(x) any(strcmpi(x,valid_values));
			
			p.addParameter('MainStruct',default_values,check_values);
			p.addParameter('ScatStruct',default_values,check_values);
			
			parse(p,varargin{:})
			
			obj.main_struct = p.Results.MainStruct;
			obj.scat_struct = p.Results.ScatStruct;
			
		end
		
		function gen_cent(obj,rows_cyl,cyl_period,arg3,varargin)
			
			if strcmpi(obj.main_struct,'circle')
				
				n_col = arg3;
				
				obj.gen_cent_circ(rows_cyl,cyl_period,n_col,varargin{:})
				
			elseif strcmpi(obj.main_struct,'rectangle')
				
				area_width = arg3;
				
				obj.gen_cent_rect(rows_cyl,cyl_period,area_width,varargin{:})
				
			end
			
		end
		
		function gen_cent_rect(obj,rows_cyl,cyl_period,area_width,varargin)
			
			p = inputParser;
			
			default_values = 'line';
			valid_values = {'line','hexagonal'};
			check_values = @(x) any(strcmpi(x,valid_values));
			
			p.addParameter('StructShape',default_values,check_values);
			
			parse(p,varargin{:})
			
			insert = @(num, arr, pos) cat(2, arr(1:pos-1), num, arr(pos:end));
			
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
			
			if strcmpi(p.Results.StructShape,'line')
			
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
		
		
		function gen_cent_circ(obj,rows_cyl,cyl_period,n_col,varargin)
			
			p = inputParser;
			
			default_values = {'-' , 'rectangle'};
			valid_values = {{'plus','+','minus','-'},{'rectangle','hexagonal'}};
			check_values = {@(x) any(strcmpi(x,valid_values{1})) , @(x) any(strcmpi(x,valid_values{2}))};
			
			p.addParameter('Scat_pm',default_values{1},check_values{1});
			p.addParameter('StructShape',default_values{2},check_values{2});
			
			parse(p,varargin{:})
			
			insert = @(num, arr, pos) cat(2, arr(1:pos-1), num, arr(pos:end));
			
			if any(strcmpi(p.Results.Scat_pm,{'plus','+'}))
				
				pm = 1;
				
			else
				
				pm = -1;
				
			end
			
			if mod(rows_cyl,2) == 0				% Checks if rows_cyl is an even number.

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
			
			if strcmpi(p.Results.StructShape,'rectangle')
			
				cyl_cent_x = zeros(1,rows_cyl);
				
			else
				
				switch mod(rows_cyl,4)
					
					case 1
						
						tot_cyl = n_col + floor(rows_cyl/4) * 4 * n_col * (1 + pm/(2*n_col));
						
					case 2
						
						tot_cyl = rows_cyl * n_col * (1 + pm/(2*n_col));
						
					case 3
						
						tot_cyl = n_col + 2*(n_col + pm) + floor(rows_cyl/4) * 4 * n_col * (1 + pm/(2*n_col));
						
					case 0
						
						tot_cyl = rows_cyl * n_col * (1 + pm/(2*n_col));
						
				end
			
				cyl_cent_x = zeros(1,tot_cyl);

				count = 1;
				
				for i = 1:rows_cyl
					
					if any(mod(i,4) == [0 , 1])
					
						cyl_cent_x(count:count+n_col-1) = cyl_period * ((-floor(n_col/2)):(floor(n_col/2)));
						
						cyl_cent_y = insert(repmat(cyl_cent_y(count),1,n_col-1),cyl_cent_y,count+1);
						
						count = count + n_col;
						
					elseif any(mod(i,4) == [2 , 3])
						
						cyl_cent_x(count:count+n_col+pm-1) = cyl_period * (-(n_col+pm-1)/2:(n_col+pm-1)/2);
						
						cyl_cent_y = insert(repmat(cyl_cent_y(count),1,n_col+pm-1),cyl_cent_y,count+1);
						
						count = count + n_col + pm;
						
					end
					
				end

				%{
				for i = 1:rows_cyl

					if any(mod(i,4) == [0 , 1])

						cyl_cent_x(count) = 0;

						count = count+1;

					elseif any(mod(i,4) == [2 , 3])

						cyl_cent_x([count,count+n_pm]) = [area_width/2 , -area_width/2];

						cyl_cent_y = insert(cyl_cent_y(count),cyl_cent_y,count+1);

						count = count+2;

					else

						error("Mod stuff doesn't work")

					end

				end
				%}
				
			end
			
			obj.cent_x = cyl_cent_x;
			obj.cent_y = cyl_cent_y;
			
		end
		
	end
	
end

