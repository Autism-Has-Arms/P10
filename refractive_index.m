
d = 100;
lambda = 700;
k_0 = 2*pi/lambda;
p = 200;
q = 400;
nr = linspace(0,6,p);
ni = linspace(0,6,q);
n = nr + 1i*ni.';
a = 1i*k_0*d*n;
r = 0.13;
t = 0.4;
f = (((exp(2*a) - 1) + sqrt( (1-exp(2*a)).^2 + 4*(r^2)*exp(2*a) ))./(2*r*exp(2*a))) - sqrt(( exp(a) - t )./( exp(a) - t*exp(2*a) ));

for l=1:p-1
    for k=1:q-1
        counter = 0;
        d1 = abs( atan(imag(f(k,l))/real(f(k,l))) - atan(imag(f(k+1,l))/real(f(k+1,l))) );
        d2 = abs( atan(imag(f(k,l))/real(f(k,l))) - atan(imag(f(k,l+1))/real(f(k,l+1))) );
        d3 = abs( atan(imag(f(k+1,l))/real(f(k+1,l))) - atan(imag(f(k+1,l+1))/real(f(k+1,l+1))) );
        d4 = abs( atan(imag(f(k,l+1))/real(f(k,l+1))) - atan(imag(f(k+1,l+1))/real(f(k+1,l+1))) );
        f(k,l);
%         d1
%         d2
%         d3
%         d4
        if d1 > pi/2
            counter = counter+1;
        end
        if d2 > pi/2
            counter = counter+1;
        end
        if d3 > pi/2
            counter = counter+1;
        end
        if d4 > pi/2
            counter = counter+1;
        end

        if counter == 1
          nr(k)
          ni(l)
          k
          l
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
		
		error('None of the conditionals selected.')
		
	end

end

