% clear all
% cac = data();

load('Struct_data.mat');
% wavelength = cac{1};
% R = cac{2};
% T = cac{3};
R = reflectance;
T = transmittance;

r_index = zeros(length(wavelength),2);

for h=1:length(wavelength)
	
	nr_int = [];
	ni_int = [];
	
    d = distance(h);
%     d = 100;
    lambda = wavelength(h);
    k0 = 2*pi/lambda;
    
    % Defines the grid for which the function is calculated.
    p = 2000;
    q = 2000;
    nr = linspace(0.001,5,p);
    ni = linspace(0.001,5,q);
    
    % n is the matrix contaning the points in the complex plane for which
    % f is calculated.
    n = nr + 1i*ni.';
    
    a = 1i*k0*d*n;
    
    % R and T are vectors containing the reflectivity and transmitivity for
    % different measurements.
     r = R(h);
     t = T(h);

    % Define the function you want to find zeros for. 
%     f = (((exp(2*a) - 1) + sqrt( (1-exp(2*a)).^2 + 4*(r^2)*exp(2*a) ))./(2*r*exp(2*a))) - sqrt(( exp(a) - t )./( exp(a) - t*exp(2*a) ));
    n1 = 1;
	r12 = (n - n1)./(n + n1);
	f = r - (r12 .* (1 - exp(2i.*k0*n*d)))./(1 - r12.^2 .* exp(2i*k0*n*d));
    
    % Phase_plot is only used to to get a view of the phase of f in the
    % comple plane nr + ni. This is not needed for the script to find the
    % zero-points.
    
    phase_plot = zeros(q,p);
	
    for l=1:p-1
        for k=1:q-1

            counter = 0;
            % The 4 d's gives the change of phase  of function f between neighbouring
            % points in the complex plane. Together, these gives the phase
            % change along each edge of a rectangle formed by 4
            % neighbouring points.
            d1 = abs( phase(f(k,l)) - phase(f(k+1,l)) );
            d2 = abs( phase(f(k,l)) - phase(f(k,l+1)) );
            d3 = abs( phase(f(k+1,l)) - phase(f(k+1,l+1)) );
            d4 = abs( phase(f(k,l+1)) - phase(f(k+1,l+1)) );
            
            % Counter starts as 0, and is increased by 1 for each edge of
            % the tringle that has a large change in phase.
            if d1 > 1.0*pi
                counter = counter+1;
            end
            if d2 > 1.0*pi
                counter = counter+1;
            end
            if d3 > 1.0*pi
                counter = counter+1;
            end
            if d4 > 1.0*pi
                counter = counter+1;
            end
            
            % If only one edge has a large change in phase, the rectangle
            % contains a singularity.
            if counter == 1 %&& abs(f(k,l)) < 0.3 && abs(nr(l)) > 0.1

				nr_int = linspace(l,l+1,p);
				ni_int = linspace(k,k+1,q);
				
				for loop_nr = 1:length(nr_int)-1
					
					for loop_ni = 1:length(ni_int)-1
						
						counter2 = 0;
						
						u1 = abs( phase(f(loop_ni,loop_nr)) - phase(f(loop_ni+1,loop_nr)) );
						u2 = abs( phase(f(loop_ni,loop_nr)) - phase(f(loop_ni,loop_nr+1)) );
						u3 = abs( phase(f(loop_ni+1,loop_nr)) - phase(f(loop_ni+1,loop_nr+1)) );
						u4 = abs( phase(f(loop_ni,loop_nr+1)) - phase(f(loop_ni+1,loop_nr+1)) );
						
						if u1 > 1.0*pi
							counter2 = counter2+1;
						end
						if u2 > 1.0*pi
							counter2 = counter2+1;
						end
						if u3 > 1.0*pi
							counter2 = counter2+1;
						end
						if u4 > 1.0*pi
							counter2 = counter2+1;
						end
						
						if counter2 == 1

							%nr(l) and ni(k) gives the real and imaginary value for n at
							%the point which is the first corner of the trianle containing
							%the singularity.


							% f of n is printed to show that the point which has been found
							% is indeed a value for which f is zero. f(k,l) should
							% therefore be a very low value.
							% f(k,l)
							% Once a point has been found is is added to a list.

							n_cal = nr(loop_nr) + 1i*ni(loop_ni)
							
						end
						
					end
					
				end
				
				if size(nr_int,1) > 1 %r_index(h) ~= 0
				
					error('Too many zero-points found')
				
				end
				
            end
            
            % Gives a matrix containing the phase of f for each nr + ni
            % value. Used for a phase_plot if this is desired, but not used
            % for finding the zero-points.
            phase_plot(k,l) = phase(f(k,l));
			
		end
		
	end
	
	r12_cal = (n_cal-n1)/(n_cal + n1);
	trans = ((1-r12_cal^2)*exp(1i*k0*d*n_cal))/(1 - r12_cal^2*exp(2i*k0*d*n_cal));
	accuracy = abs(trans - t)/abs(t);
	r_index(h,1) = n_cal;
	r_index(h,2) = accuracy;
	
end

% surf(nr,ni,phase_plot)
% shading interp
% view(2)
% calculated the phase for a given complex number. Note that the phase goes
% from -pi/2 to 3pi/2. For a complex number = 0, the phase is undefined,
% and the functions returns an error.
function f = phase(f)

	if real(f) == 0 && imag(f) > 0
		
		f = pi/2;
		
	elseif real(f) == 0 && imag(f) < 0
		
		f = -pi/2;
		
	elseif real(f) > 0 && imag(f) < 0
		
		f = atan(imag(f)/real(f)) + 2*pi;
    
    elseif real(f) > 0 && imag(f) > 0
        
        f = atan(imag(f)/real(f));
		
	elseif real(f) < 0
		
		f = atan(imag(f)/real(f)) + pi;
		
    else
		

		error('None of the conditionals selected.')
		
	end

end


function    cac = data()
        str = fileread( 'testdata1.txt' );
        str = strrep( str, ' ', '' );
        cac = textscan( str, '%f%f%f%f', 'Delimiter', ',' ); 
end
