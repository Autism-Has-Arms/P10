classdef zone_determination < handle
	%
	
	properties
		
		PML
		env
		upper
		lower
		scat
		scat_whole
		scat_split
		
	end
	
	
	methods
		
		function zone_det(obj,dl,bt,varargin)
			
			p = inputParser;
			
			default_values = 'circle';
			valid_values = {'circle','rectangle'};
			check_values = @(x) any(strcmpi(x,valid_values));

			addParameter(p,'ScatType',default_values,check_values);
			
			parse(p,varargin{:})
			
			%% Environment faces
			
			obj.env = find(sum(bt(:,1:(end - 1)),2) - 1 == 1).';
			
			if length(obj.env) ~= (2)
				
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
			
			%% Scatterer faces
			
			scat_geom = bt(:,(3):(end - 1));
			
			if exist('key','var')
				
				scat_geom(obj.env,:) = 0;
				
			end
			
			[row,col] = ind2sub(size(scat_geom),find(scat_geom));
			ind_cell = mat2cell(col == unique(col).',length(col),ones(1,length(unique(col))));
			
			for i = 1:length(ind_cell)
				
				scat = row(ind_cell{i});
				
				if length(scat) == 1
					
					obj.scat_whole = [obj.scat_whole scat];
					
				else
					
					if bt(scat(1),end)
						
						scat = flipud(scat);
						
					end
					
					obj.scat_split = [obj.scat_split scat];
					
				end
				
			end
			
			obj.scat = [obj.scat_whole reshape(obj.scat_split.',1,[])];
			
			%% Perfectly Matched Layer
			
			obj.PML = find(sum(bt(:,1:end - 1),2) == 1).';

			if bt(obj.PML(1),end)

				obj.PML = fliplr(obj.PML);

			end
				
			obj.upper = obj.env(bt(obj.env,end) == 0);
			obj.lower = obj.env(bt(obj.env,end) == 1);

			obj.upper = [obj.upper , obj.PML(bt(obj.PML,end) == 0)];
			obj.lower = [obj.lower , obj.PML(bt(obj.PML,end) == 1)];
			
		end
		
	end
	
end

