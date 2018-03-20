function pct = disppct(varargin)

    if nargin == 0
        
        dispstat('','init');
        dispstat(sprintf('Beginning the process.'),'keepthis','timestamp');
        pct = -1;
        
    elseif nargin == 4

        count = varargin{1};
        lngt = varargin{2};
        pct = varargin{3};
        n = varargin{4};
        
        pct1 = floor(count/lngt.*100);

        if pct ~= pct1

            dispstat(['For-loop #' num2str(n) ': ' num2str(pct1) '%'],'timestamp');

            pct = pct1;

        end

        if count == lngt
            
%             if n == 2
% 
%                 dispstat('Finished.','keepprev','timestamp');
%                 
%             else
                
                dispstat(' ','keepprev');
				
				pct = 0;
                
%             end

        end

    end

end