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
    p = 500;
    q = 500;
    nr = linspace(0.1,6,p);
    ni = linspace(0.1,6,q);
    n = nr + 1i*ni.';
    a = 1i*k_0*d*n;
    r = R(h);
    t = T(h);
    f = (((exp(2*a) - 1) + sqrt( (1-exp(2*a)).^2 + 4*(r^2)*exp(2*a) ))./(2*r*exp(2*a))) - sqrt(( exp(a) - t )./( exp(a) - t*exp(2*a) ));
    phase_plot = zeros(q,p);
    for l=1:p-1
        for k=1:q-1

            counter = 0;
            d1 = abs( phase(f(k,l)) - phase(f(k+1,l)) );
            d2 = abs( phase(f(k,l)) - phase(f(k,l+1)) );
            d3 = abs( phase(f(k+1,l)) - phase(f(k+1,l+1)) );
            d4 = abs( phase(f(k,l+1)) - phase(f(k+1,l+1)) );
            f(k,l);
%         d1
%         d2
%         d3
%         d4
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

            if counter == 1
            h
            k
            l
            nr(l)
            ni(k)
            f(k,l)
            r_index(h,1) = r_index(h,1) + nr(l) + 1i*ni(k);
            end
            phase_plot(k,l) = phase(f(k,l));
        end
    end
end

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
		
%         f = pi/2;
		error('None of the conditionals selected.')
		
	end

end


function    cac = data()
        str = fileread( 'testdata1.txt' );
        str = strrep( str, ' ', '' );
        cac = textscan( str, '%f%f%f%f', 'Delimiter', ',' ); 
end
