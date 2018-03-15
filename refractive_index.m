cac = data();

% data = load('testdata1.txt');
wavelength = cac{1};
R = cac{2};
T = cac{3};

r_index = zeros(length(wavelength),1);

for h=3:4
    d = 100;
    lambda = wavelength(h);
    k_0 = 2*pi/lambda;
    
    % Defines the grid for which the function is calculated.
    p = 500;
    q = 500;
    nr = linspace(0.1,6,p);
    ni = linspace(0.1,6,q);
    
    % n is the matrix contaning the points in the complex plane for which
    % f is calculated.
    n = nr + 1i*ni.';
    
    a = 1i*k_0*d*n;
    
    % R and T are vectors containing the reflectivity and transmitivity for
    % different measurements.
    r = R(h);
    t = T(h);
    
    % Define the function you want to find zeros for. 
    f = (((exp(2*a) - 1) + sqrt( (1-exp(2*a)).^2 + 4*(r^2)*exp(2*a) ))./(2*r*exp(2*a))) - sqrt(( exp(a) - t )./( exp(a) - t*exp(2*a) ));
    
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
            if d1 > pi
                counter = counter+1;
            end
            if d2 > pi
                counter = counter+1;
            end
            if d3 > pi
                counter = counter+1;
            end
            if d4 > pi
                counter = counter+1;
            end
            
            % If only one edge has a large change in phase, the rectangle
            % contains a singularity.
            if counter == 1
            h
            k
            l
            
            %nr(k) and ni(l) gives the real and imaginary value for n at
            %the point which is the first corner of the trianle containing
            %the singularity.
            nr(k)
            ni(l)
            
            % f of n is printed to show that the point which has been found
            % is indeed a value for which f is zero. f(k,l) should
            % therefore be a very low value.
            f(k,l)
            
            % Once a point has been found is is added to a list.
            r_index(h,1) = r_index(h,1) + nr(l) + 1i*ni(k);
            end
            
            % Gives a matrix containing the phase of f for each nr + ni
            % value. Used for a phase_plot if this is desired, but not used
            % for finding the zero-points.
            phase_plot(k,l) = phase(f(k,l));
        end
    end
end

% calculated the phase for a given complex number. Note that the phase goes
% from -pi/2 to 3pi/2. For a complex number = 0, the phase is undefined,
% and the functions returns an error.
function f = phase(f)

	if real(f) == 0 && imag(f) > 0
		
		f = pi/2;
		
	elseif real(f) == 0 && imag(f) < 0
		
		f = -pi/2;
		
	elseif real(f) > 0
		
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