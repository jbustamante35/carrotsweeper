function [wave] = buildWave(wave_tensor,v_vec,sig_vec,theta_val,rho_val,R,NF,NP)
    oPath = '/home/nate/labTMP/';
    mkdir(oPath);
    wave = zeros(size(theta_val));
   
    
    dr = 1;
    buildR = 0:max(rho_val);
    K = 0:(NF-1);
    [WG(:,:,1),WG(:,:,2)] = ndgrid(linspace(R(1),R(2),NP(1)),K);
    
    %U = zeros(size(ORG));
    
    v = VideoWriter([oPath 'wave1.mpeg'],'Grayscale AVI');
    open(v)
    
    for r = buildR
        fidx = find(rho_val >= r & rho_val < r + dr);
        
        %rhoVec = rho_val(fidx);
        %K = 0:10;
        
        %K = 0;
        rhoVec = r;
        for k = K
            % for disk square -debug
            %theta_value = linspace(0,2*pi,size(ORG,2));
            
            
            [iT] = iWaveTensor(WG,rhoVec,k,wave_tensor);
            
            theta_value = theta_val(fidx) + pi;
            %iT(:,2) = 0;
            f = makeWave(iT,k,theta_value,v_vec(k+1),sig_vec(k+1));
            wave(fidx) = wave(fidx) + f;
            %f = makeWave(iT,k,theta_value,size(ORG,2));
            %U(r+1,:) = U(r+1,:)+f;
            
        end
        imshow(real(wave),[])
        writeVideo(v,bindVec(real(wave)));
        
        %drawnow
    end
    
    close(v)
    %{
    for rho = 1:size(wave_tensor,1)
        
        
        for w = 1:size(wave_tensor,2)
            
            wave = wave + wave_tensor(rho,w,1)*exp(1i*omega(w-1,wave_tensor(rho,w,2),v_vec(rho,w),sig_vec(rho,w),theta_val));
        
        end
        
        
        imshow(wave,[]);
        drawnow
    end
    %}
end

function [iT] = iWaveTensor(WG,rho,k,wave_tensor)
    amp = interp2(WG(:,:,2),WG(:,:,1),wave_tensor(:,:,1),k*ones(size(rho)),rho,'linear');
    ang = interp2(WG(:,:,2),WG(:,:,1),wave_tensor(:,:,2),k*ones(size(rho)),rho,'linear');
    iT = cat(2,amp,ang);
end


function [wave] = makeWave(wave_tensor,kValue,thetaValue,VEL,SIG)
    %omega = @(k,phi,v,sig,theta)2*pi*N^-1*k*(theta + v*sig) + phi;
    omega = @(k,phi,v,sig,theta)k*(theta + v*sig) + phi;
    wave = wave_tensor(:,1).*exp(1i*omega(kValue,wave_tensor(:,2),VEL,SIG,thetaValue));
end