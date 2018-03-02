
d = 100;
lambda = 700;
k_0 = 2*pi/lambda;
p = 2000;
q = 10000;
nr = linspace(0,6,p);
ni = linspace(0,6,2*q);
n = nr + 1i*ni';
a = 1i*k_0*d*n;
r = 0.13;
t = 0.4;
f = ((exp(2*a) - 1) + sqrt( (1-exp(2*a)).^2 + 4*r^2*exp(2*1) ))/(2*r*exp(2*a)) - sqrt(( exp(a) - t )/( exp(a) - t*exp(2*a) ));

for k=1:p-1
    for l=1:q-1
        counter = 0;
        d1 = abs( atan(imag(f(k,l))/real(f(k,l))) - atan(imag(f(k+1,l))/real(f(k+1,l))) );
        d2 = abs( atan(imag(f(k,l))/real(f(k,l))) - atan(imag(f(k,l+1))/real(f(k,l+1))) );
        d3 = abs( atan(imag(f(k+1,l))/real(f(k+1,l))) - atan(imag(f(k+1,l+1))/real(f(k+1,l+1))) );
        d4 = abs( atan(imag(f(k,l+1))/real(f(k,l+1))) - atan(imag(f(k+1,l+1))/real(f(k+1,l+1))) );
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

% for k=1:q
%     for l=1:q
%         if f(k,l)<0.01
%             nr(k)
%             ni(l)
%         end
%     end
% end

% surf(nr,ni,f)
% shading interp
