function [f] = lineSample(stack,sheets)
    f = zeros(size(sheets,1),size(sheets,2),size(stack,3),size(stack,4));
    for fr = 1:size(stack,4)
        for k = 1:size(stack,3)
            f(:,:,k,fr) = interp2(stack(:,:,k,fr),sheets(:,:,1),sheets(:,:,2));
        end
    end
end


%{
        sheetF = squeeze(f(:,s,k,:))';
        E = edge(sheetF);
        [g1 g2] = gradient(sheetF);
        %vel = g2(find(E)).*g1(find(E)).^-1;
        A = atan2(g1(find(E)),g2(find(E)));
        A(A<-pi/2) = A(A<-pi/2) + pi;
        A(A>pi/2) = A(A>pi/2) - pi;
        N1 = [g1(find(E));-g1(find(E))];
        N2 = [g2(find(E));-g2(find(E))];
        [f,xi] = ksdensity(A);
        figure(h3)
        plot(f);
        [~,midx]= max(f);
        g(s) = mean(A);
        g(s) = xi(midx);
        %{
        figure(h3)
       
        %plot(N1,N2,'.')
        %drawnow
        %g(s) = mean().^2 + g2(find(E)).^2);
      
        fprintf(['done with line sample:' num2str(s) ':' num2str(size(sheets,2)) '\n']);
        figure(h1)
        imshow(f(:,:,:,s),[0 1]);
        figure(h2)
        imshow(stack(:,:,:,1),[0 1])
        hold on
        plot(sheets(:,s,1),sheets(:,s,2),'r');
        hold off
        drawnow
        %}
    %end
    [~,midx] = max(g);
end
%}