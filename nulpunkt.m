Solution to problem 2:
%Calculate mode indices of asymmetric dielectric waveguides using bisection
 
clear all; close all
n1=1.0;
n2=3.5;
n3=1.45;
 
%Waveguide_thicknesses
d_start = 0; 
d_end   = 1000; 
N_d = 200; %Number of waveguide thicknesses
 
wavelength=1500; %[nm]
 
s_polarization = true;
p_polarization = false;
 
sampling_points = 1000;
Number_of_mode_indices_to_search_for = 5;
accuracy = 1e-6;
 
plot_f = false;
plot_phase = false;
plot_nm_vs_d = true;
plot_field = false;
 
%----- program -------
output=zeros(N_d,Number_of_mode_indices_to_search_for);
eps1=n1^2; eps2=n2^2; eps3=n3^2;
k0=2*pi/wavelength
 
ky1=@(nm) sqrt(-k0^2*nm.^2+k0^2*eps1);
ky2=@(nm) sqrt(-k0^2*nm.^2+k0^2*eps2);
ky3=@(nm) sqrt(-k0^2*nm.^2+k0^2*eps3);
if p_polarization == true,
   r12=@(nm) (eps2*ky1(nm)-eps1*ky2(nm))./(eps2*ky1(nm)+eps1*ky2(nm));
   r23=@(nm) (eps3*ky2(nm)-eps2*ky3(nm))./(eps3*ky2(nm)+eps2*ky3(nm));
end
if s_polarization == true,
   r12=@(nm) (ky1(nm)-ky2(nm))./(ky1(nm)+ky2(nm));
   r23=@(nm) (ky2(nm)-ky3(nm))./(ky2(nm)+ky3(nm));
end
 
dv=linspace(d_start,d_end,N_d);
nm_result=0;
for j_d=1:N_d
    d=dv(j_d);
    if n1==n3,
        f=@(nm) ((r12(nm).^(-1)).*exp(-1i*ky2(nm)*d)-...
            r12(nm).*exp(1i*ky2(nm)*d));
    else
        f=@(nm) (1+r12(nm).*r23(nm).*exp(2i*ky2(nm)*d));
    end
    
    nmax=n2; nmin=max([n1 n3]);
    nmv=linspace(nmin+1e-6,nmax-1e-6,sampling_points);
    
    if plot_f==true,
        plot(nmv,real(f(nmv)),'k','linewidth',2)
        hold on
        plot(nmv,imag(f(nmv)),'k--','linewidth',2)
        hold off
        set(gca,'fontsize',28,'FontName','Times New Roman')
        xlabel('Mode index','fontsize',28,'FontName','Times New Roman')
        h=legend('Real\{f({\it{n}}_m)\}','Imag\{f({\it{n}}_m)\}'); 
        set(h,'box','off')
        set(h,'FontName','Times New Roman','fontsize',28)
        text(1.5,-1.3,'{\it{d}} = 150 nm','fontsize',28,...
            'FontName','Times New Roman') 
        text(1.5,-0.8,'Air-Silicon-Quartz','fontsize',28,...
            'FontName','Times New Roman')
        axis([nmin nmax -2.1 2.1])
        set(gca,'linewidth',1)
    end
    if plot_phase==true,
        plot(nmv,phase(f(nmv)),'k','linewidth',2)
        set(gca,'fontsize',28,'FontName','Times New Roman')
        xlabel('Mode index','fontsize',28,'FontName','Times New Roman')
        ylabel('Phase [rad]','fontsize',28,'FontName','Times New Roman')
        text(1.5,0.85,'{\it{d}} = 150 nm','fontsize',28,...
            'FontName','Times New Roman')
        text(1.5,1.2,'Air-Silicon-Quartz','fontsize',28,...
            'FontName','Times New Roman')
        axis([nmin nmax -pi/2-0.1 pi/2+0.1])
        set(gca,'linewidth',1)
    end
    
    %Find mode index using bisection
    mode_index_counter = 1;
    for j=1:sampling_points-1,
       if abs(phase(f(nmv(j)))-phase(f(nmv(j+1))))>pi*3/4,
          %There is a solution in the interval
          nma=nmv(j); nmb=nmv(j+1);
          fa=f(nma); fb=f(nmb);
          while nmb-nma>accuracy,
             nmc=(nma+nmb)/2; fc=f(nmc);
             if abs(phase(fa)-phase(fc))>pi*3/4, 
                 nmb=nmc; fb=fc;
             else
                 nma=nmc; fa=fc;
             end
          end
          nm_result(mode_index_counter)=(nma+nmb)/2;
          mode_index_counter=mode_index_counter+1;
       end
    end
    nm_result
    dummy=sort(nm_result,'descend');
    if length(dummy)>=Number_of_mode_indices_to_search_for,
        output(j_d,:)=dummy(1:Number_of_mode_indices_to_search_for);
    else
        output(j_d,1:length(dummy))=dummy;
    end
end
 
if plot_nm_vs_d == true,
    for j=1:Number_of_mode_indices_to_search_for,
        index1=find(output(:,j)>0,1);
        index2=N_d;
        plot(dv(index1:index2),output(index1:index2,j),'k','linewidth',2)
        hold on
    end
    axis([d_start d_end n3 n2])
    set(gca,'fontsize',28,'FontName','Times New Roman');
    xlabel('{\it{d}} [nm]','fontsize',28,'FontName','Times New Roman')
    ylabel('Mode index','fontsize',28,'FontName','Times New Roman')
    set(gca,'linewidth',1,'FontName','Times New Roman')
end
 
if plot_field == true,
   A=1;
   y1v=linspace(0,1000,1000);
   kx=k0*output(1,1);
   ky1=sqrt(k0^2*eps1-kx^2);
   ky2=sqrt(k0^2*eps2-kx^2);
   ky3=sqrt(k0^2*eps3-kx^2);
   E1v=exp(1i*ky1*y1v);
   if s_polarization == true,
      A=0.5*(1+ky1/ky2); B=0.5*(1-ky1/ky2); 
   end
   if p_polarization == true,
      A=0.5*(1+ky1/ky2*eps2/eps1); B=0.5*(1-ky1/ky2*eps2/eps1);  
   end
   C=A*exp(1i*ky2*(-d))+B*exp(-1i*ky2*(-d));
   
   y2v=linspace(-d,0,200);
   E2v=A*exp(1i*ky2*y2v)+B*exp(-1i*ky2*y2v);
   
   y3v=linspace(-1000-d,-d,1000);
   E3v=C*exp(-1i*ky3*(d+y3v));
   
   yv=[y3v y2v y1v]; Ev=[E3v E2v E1v];
   plot(yv,real(Ev)/max(real(Ev)),'k','linewidth',2)
   set(gca,'fontsize',28,'FontName','Times New Roman')
   xlabel('{\it{y}} [nm]','fontsize',28,'FontName','Times New Roman')
   ylabel('Electric field','fontsize',28,'FontName','Times New Roman')
   axis([-1000 500 -1.1 1.1])
   
   hold on
   A=1;
   y1v=linspace(0,1000,1000);
   kx=k0*output(1,2);
   ky1=sqrt(k0^2*eps1-kx^2);
   ky2=sqrt(k0^2*eps2-kx^2);
   ky3=sqrt(k0^2*eps3-kx^2);
   E1v=exp(1i*ky1*y1v);
   if s_polarization == true,
      A=0.5*(1+ky1/ky2); B=0.5*(1-ky1/ky2); 
   end
   if p_polarization == true,
      A=0.5*(1+ky1/ky2*eps2/eps1); B=0.5*(1-ky1/ky2*eps2/eps1);  
   end
   C=A*exp(1i*ky2*(-d))+B*exp(-1i*ky2*(-d));
   
   y2v=linspace(-d,0,200);
   E2v=A*exp(1i*ky2*y2v)+B*exp(-1i*ky2*y2v);
   y3v=linspace(-1000-d,-d,1000);
   E3v=C*exp(-1i*ky3*(d+y3v));
   yv=[y3v y2v y1v]; Ev=[E3v E2v E1v];
   plot(yv,real(Ev)/max(real(Ev)),'k--','linewidth',2)
   
   hold on
   h=legend('Mode 1','Mode 2'); set(h,'FontName','Times New Roman');
   set(h,'box','off')
   x=[-d -d]; y=[-1.1 1.1];
   plot(x,y,'k:','linewidth',2)
   x=[0 0]; y=[-1.1 1.1];
   plot(x,y,'k:','linewidth',2)
   
   text(100,-0.5,'{\it{n}} = 1','fontsize',28,'FontName','Times New Roman')
   text(-390,-0.5,'{\it{n}} = 3.5','fontsize',28,...
       'FontName','Times New Roman')
   text(-900,-0.5,'{\it{n}} = 1.45','fontsize',28,...
       'FontName','Times New Roman')
   set(gca,'linewidth',1)
end

Solution to problem 3:
%Calculates the modes of a gold-air-gold waveguide using Newton-Raphson
 
clear all; close all
 
wavelength=800; %[nm]
n1=0.1532 + 4.8984i; %refractive index of gold at wavelength 800 nm
n2=1.0;
n3=n1; 
 
%Waveguide_thicknesses
d_start = 0.3; %[nm]
d_end   = 100; %[nm]
N_d = 200; %Number of waveguide thicknesses
 
p_polarization = true;
 
nr_min=1.01; 
nr_max=40; 
ni_min=1e-7;
ni_max=2.5; 
 
sampling_points_nr = 200; %number of intervals along real axis
sampling_points_ni = 200; %number of intervals along imaginary axis
 
accuracy = 1e-6;
Number_of_modes_to_search_for = 5;
 
 
plot_f = false;
plot_phase = false;
plot_nm_vs_d = true;
plot_field = false;
 
 
%----- program -------
output=zeros(N_d,Number_of_modes_to_search_for);
eps1=n1^2; eps2=n2^2; eps3=n3^2;
k0=2*pi/wavelength
 
%We require that the imaginary part of $k_{y,i}$ >= 0
ky1d=@(nm) sqrt(-k0^2*nm.^2+k0^2*eps1)
ky1=@(nm) ky1d(nm).*(imag(ky1d(nm))>=0) - ky1d(nm).*(imag(ky1d(nm))<0);
ky2d=@(nm) sqrt(-k0^2*nm.^2+k0^2*eps2)
ky2=@(nm) ky2d(nm).*(imag(ky2d(nm))>=0) - ky2d(nm).*(imag(ky2d(nm))<0);
ky3d=@(nm) sqrt(-k0^2*nm.^2+k0^2*eps3)
ky3=@(nm) ky3d(nm).*(imag(ky3d(nm))>=0) - ky3d(nm).*(imag(ky3d(nm))<0);
r12=@(nm) (eps2*ky1(nm)-eps1*ky2(nm))./(eps2*ky1(nm)+eps1*ky2(nm));
r23=@(nm) (eps3*ky2(nm)-eps2*ky3(nm))./(eps3*ky2(nm)+eps2*ky3(nm));
r12m=@(nm) 2*k0*nm*k0^2*eps1*eps2*(eps1-eps2)./ky1(nm)./ky2(nm)./((eps2*ky1(nm)+eps1*ky2(nm)).^2);
r23m=@(nm) 2*k0*nm*k0^2*eps2*eps3*(eps2-eps3)./ky2(nm)./ky3(nm)./((eps3*ky2(nm)+eps2*ky3(nm)).^2);
 
dv=linspace(d_start,d_end,N_d);
for j_d=1:N_d
    d=dv(j_d)
    %We define the function f and f'
    f=@(nm) (1+r12(nm).*r23(nm).*exp(2i*ky2(nm)*d));
    fm=@(nm) (r12m(nm).*r23(nm)+r12(nm).*r23m(nm)+...
        r12(nm).*r23(nm).*2i.*(-k0.*d.*nm./ky2(nm))).*exp(2i*ky2(nm)*d)*k0;
    
    %Calculates f and the phase of f on a grid of points
    nmrv=linspace(nr_min,nr_max,sampling_points_nr);
    nmiv=linspace(ni_min,ni_max,sampling_points_ni);
    for jnr=1:sampling_points_nr,
        for jni=1:sampling_points_ni,
            nm=nmrv(jnr)+1i*nmiv(jni);
            fM(jni,jnr)=f(nm);
            phaseM(jni,jnr)=phase(fM(jni,jnr));
            ky1M(jni,jnr)=ky1(nm);
            ky2M(jni,jnr)=ky2(nm);            
            ky3M(jni,jnr)=ky3(nm);
        end
    end
    
    if plot_f==true,
        contour(nmrv,nmiv,real(fM),[0 0],'k','linewidth',2)
        hold on
        contour(nmrv,nmiv,imag(fM),[0 0],'k--','linewidth',2)
        
        contour(nmrv,nmiv,real(ky1M),[0 0],'k:','linewidth',1)
        contour(nmrv,nmiv,real(ky2M),[0 0],'k:','linewidth',1)
        contour(nmrv,nmiv,real(ky3M),[0 0],'k:','linewidth',1)        
        
        hold off
        set(gca,'fontsize',28,'FontName','Times New Roman')
        xlabel('Real part of mode index','fontsize',28,...
            'FontName','Times New Roman')
        ylabel('Imaginary part of mode index','fontsize',28,...
            'FontName','Times New Roman')
        h=legend('Real\{f({\it{n}}_m)\} = 0','Imag\{f({\it{n}}_m)\} = 0'); 
        set(h,'box','off')
        set(h,'FontName','Times New Roman','fontsize',28)
        set(gca,'linewidth',1)
        axis image
    end
    if plot_phase==true,
        pcolor(nmrv,nmiv,phaseM), shading interp, h=colorbar;
        axis image, colormap(gray)
        hold off
        set(gca,'fontsize',28,'FontName','Times New Roman')
        xlabel('Real part of mode index','fontsize',28,...
            'FontName','Times New Roman')
        ylabel('Imaginary part of mode index','fontsize',28,...
            'FontName','Times New Roman')
        set(gca,'linewidth',1)
    end
    
    %Find mode index using Newton-Raphson
    nm_result=0;
    mode_index_counter = 1;
    for jnr=1:sampling_points_nr-1,
        for jni=1:sampling_points_ni-1,
            if mode_index_counter<=Number_of_modes_to_search_for,
                njumps=0;
                if abs(phaseM(jni,jnr+1)-phaseM(jni,jnr))>pi, 
                    njumps=njumps+1; end
                if abs(phaseM(jni+1,jnr+1)-phaseM(jni,jnr+1))>pi, 
                    njumps=njumps+1; end
                if abs(phaseM(jni+1,jnr)-phaseM(jni+1,jnr+1))>pi, 
                    njumps=njumps+1; end
                if abs(phaseM(jni,jnr)-phaseM(jni+1,jnr))>pi, 
                    njumps=njumps+1; end            
 
                if njumps==1, 
                    %There is a solution in the rectangle
                    nma=nmrv(jnr)+1i*nmiv(jni);
                    nmb=nmrv(jnr+1)+1i*nmiv(jni); 
                    nmc=nmrv(jnr+1)+1i*nmiv(jni+1); 
                    nmd=nmrv(jnr)+1i*nmiv(jni+1);
                    nm1=(nma+nmb+nmc+nmd)/4 %start-point for Newton-Raphson
 
                    iterations=1;
                    error=1;
                    stop=0;
                    while error>accuracy && stop==0,
                        f1=f(nm1);
                        fm1=fm(nm1);
                        nm2=nm1-f1/fm1;
                        error=abs(nm2-nm1);
                        if error>max([abs(nmb-nma)*2, abs(nmc-nmb)*2]),
                            stop=1;
                        end
                        iterations=iterations+1;
                        if iterations>100,
                            stop=1;
                        end
                        nm1=nm2;
                    end
                    if stop == 0,
                       %We found a solution, and if we have not already 
                       %found the solution one we add it to a vector
                       if max(abs(nm_result-nm1))>accuracy,
                           nm_result(mode_index_counter)=nm1; 
                           mode_index_counter=mode_index_counter+1;
                       end
                    end
                end
            end
        end
    end
    nm_result
    if length(nm_result)>=Number_of_modes_to_search_for,
        output(j_d,:)=nm_result(1:Number_of_modes_to_search_for);
    else
        output(j_d,1:length(nm_result))=nm_result;
    end
end
 
if plot_nm_vs_d == true,
    for j=1:Number_of_modes_to_search_for,
        index1=find(output(:,j)>0,1);
        index2=N_d;
        loglog(dv(index1:index2),real(output(index1:index2,j)),...
            'k','linewidth',2)
        hold on
        loglog(dv(index1:index2),imag(output(index1:index2,j)),...
            'k--','linewidth',2)
    end
    axis([d_start d_end 0.5e-2 40])
    set(gca,'fontsize',28,'FontName','Times New Roman');
    xlabel('{\it{d}} [nm]','fontsize',28,'FontName','Times New Roman')
    ylabel('Mode index','fontsize',28,'FontName','Times New Roman')
    set(gca,'linewidth',1,'FontName','Times New Roman')
    h=legend('Real part','Imaginary part')
    set(h,'fontsize',28,'FontName','Times New Roman')    
end
 
if plot_field == true,
   A=1;
   y1v=linspace(0,100,200);
   kx=k0*output(1,1);
   ky1=sqrt(k0^2*eps1-kx^2);
   ky2=sqrt(k0^2*eps2-kx^2);
   ky3=sqrt(k0^2*eps3-kx^2);
   E1v=exp(1i*ky1*y1v);
   A=0.5*(1+ky1/ky2*eps2/eps1); B=0.5*(1-ky1/ky2*eps2/eps1);  
   C=A*exp(1i*ky2*(-d))+B*exp(-1i*ky2*(-d));
   y2v=linspace(-d,0,200);
   E2v=A*exp(1i*ky2*y2v)+B*exp(-1i*ky2*y2v);
   y3v=linspace(-100-d,-d,200);
   E3v=C*exp(-1i*ky3*(d+y3v));
   yv=[y3v y2v y1v]; Ev=[E3v E2v E1v];
   plot(yv,real(Ev)/max(real(Ev)),'k','linewidth',2)
   hold on
   xv=[0 0]; yv=[0 1.1];
   plot(xv,yv,'k:','linewidth',2)
   xv=[-d -d];
   plot(xv,yv,'k:','linewidth',2)
   
   set(gca,'fontsize',28,'FontName','Times New Roman')
   xlabel('{\it{y}} [nm]','fontsize',28,'FontName','Times New Roman')
   ylabel('Magnetic field','fontsize',28,'FontName','Times New Roman')
   axis([-100-d 100 0 1.1])
   
   text(50,0.5,'Gold','fontsize',28,'FontName','Times New Roman')
   text(1,0.5,'Air','fontsize',28,'FontName','Times New Roman')
   text(-40,0.5,'Gold','fontsize',28,'FontName','Times New Roman')
   set(gca,'linewidth',1)
end
