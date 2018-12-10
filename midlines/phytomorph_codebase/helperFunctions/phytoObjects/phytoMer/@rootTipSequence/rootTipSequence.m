classdef rootTipSequence < myHS_X

    properties
        % vector field(s)        
        position_field_midline;
    end
    
    methods
        function [obj] = rootTipSequence()
            % call super
            obj = obj@myHS_X('rootTip');
            % vector field(s)
            obj.position_field_midline = phytoAflow_1();
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set function
        function [] = set(obj,fieldName,data)
            obj.(fieldName) = data;
        end
        % set data function
        function [] = setData(obj,fieldName,data)
            obj.(fieldName).setData(data);
        end
        % get function
        function [data] = get(obj,fieldName)
            data = obj.(fieldName);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % steady state
        function [profile] = steadyStateFlowAnalysis(obj,para,frameIdx)
            analysisType = 'logistic';
            switch analysisType
                case 'logistic'
                    % get points from position flow field
                    sz = obj.position_field_midline.baseSize();
                    if nargin == 2
                        pointList = obj.position_field_midline(:,:);
                        pointList = reshape(pointList,[sz size(pointList,2)]);
                    else
                        pointList = obj.position_field_midline(:,frameIdx);
                        pointList = reshape(pointList,[sz(1) numel(frameIdx) size(pointList,2)]);
                    end
                    % center the tip
                    pointList = bsxfun(@minus,pointList,pointList(1,:,:));            
                    dX = diff(pointList,1,1);
                    dV = diff(pointList,1,2);
                    dL = sum(dX.*dX,3).^.5;            
                    L = cumsum(dL,1);
                    L = L(:,1:end-1);
                    % obtain the location along the midline
                    ndX = bsxfun(@times,dX,dL.^-1);
                    dV_TM = sum(ndX(:,1:end-1,:).*dV(1:end-1,:,:),3);
                    dVT_TM = sum(dV.*dV,3).^.5;
                    dVT_TM = dVT_TM(1:end-1,:);
                    
                    
                    % create velocity model
                    %velocity_model = rootTipSequence.fitStrain(L(:),dV_TM(:),@()splineFit(4,4));
                    velocity_model = rootTipSequence.fitStrain(L(:),dV_TM(:),@()logisticFit());
                    %lDomain = linspace(min(L(:)),max(L(:)),1000);
                    lDomain = linspace(0,2000,2000);
                    [velocityProfile,strainProfile] = rootTipSequence.evalStrainModel(velocity_model,lDomain);
                    
                    
                    % create kinematics profile
                    profile.velocityProfile = velocityProfile;
                    profile.strainProfile = strainProfile;
                    profile.lDomain = lDomain;
                    profile.velocity_model = velocity_model;
                    profile.raw_domain = L(:);
                    profile.raw_velocity = dV_TM(:);
                    profile.type = 'steadyState';
                
                    % generate metrics for the strain profile
                    profile = rootTipSequence.generateStrainMetrics(profile,para);
                case 'bayes'
                    %
                    yL = repmat(1:size(L,2),[size(L,1) 1]);
                    [x1 x2] = ndgrid(40:1:1500,1:1:size(L,2));
                    vq = griddata(L,yL,dV_TM,x1,x2,'cubic');
            end
        end
        % non steady state
        function [profile] = nonSteadyStateFlowAnalysis(obj,para)
            sz = obj.position_field_midline.baseSize();
            % fit models
            for tm = 1:(sz(end)-para.width)
                profile.profileSequence(tm) = obj.steadyStateFlowAnalysis(para,tm:tm+para.width);
                %profile.profileSequence(tm) = rootTipSequence.generateStrainMetrics(profile.profileSequence(tm));
            end
            
            %{
            % eval models            
            profile.spaceDomain = [];
            profile.timeDomain = [];
            TIME = (sz(end)-width);
            for tm = 1:TIME
                SPACE = numel(profile.profileSequence{tm}.lDomain);
                profile.spaceDomain = [profile.spaceDomain;profile.profileSequence{tm}.lDomain];
                profile.timeDomain = [profile.timeDomain;tm*ones(1,SPACE)];                
                [velocityProfile(tm,:),strainProfile(tm,:)] = ...
                    rootTipSequence.evalStrainModel(profile.profileSequence{tm}.velocity_model,profile.profileSequence{tm}.lDomain);
            end
            %}
            
            % eval models
            NP = 3000;
            lDomain = linspace(0,3000,3000);            
            for tm = 1:(sz(end)-para.width)
                % eval model over domain
                [velocityProfile(tm,:),strainProfile(tm,:)] = ...
                    rootTipSequence.evalStrainModel(profile.profileSequence(tm).velocity_model,lDomain);
                % reassign velocity profile
                profile.profileSequence(tm).velocityProfile = velocityProfile(tm,:);
                profile.profileSequence(tm).strainProfile = strainProfile(tm,:);
                profile.profileSequence(tm).lDomain = lDomain;
                profile.profileSequence(tm) = rootTipSequence.generateStrainMetrics(profile.profileSequence(tm),para);
            end
            
            
            
            
            [profile.timeDomain profile.spaceDomain] = ndgrid(1:(sz(end)-para.width),1:NP);
            
            profile.velocityProfile = velocityProfile;
            profile.strainProfile = strainProfile;
            profile.type = 'nonSteadyState';
        end
        % create flow profile
        function [profile] = create_flowProfile(obj,para)
            switch para.kinematicsType
                case 'nonSteadyState'
                    profile = obj.nonSteadyStateFlowAnalysis(para);
                case 'steadyState'
                    profile = obj.steadyStateFlowAnalysis(para);
            end
        end
        % morphometric profile
        function [profile] = morphometricProfile(obj,level)
            for e = 1:numel(obj.S)
                ele = obj.S{e};
                profile.lengthProfile(e) = ele.midline_1.length();
                TM = ele.midline_1.generateTangentSpace_atP(1,level);
                tmp = TM.d;
                tmp(1:2) = -(tmp(1:2));
                tmp = reshape(tmp(1:4),[2 2]);
                tmp = phytoAaffine(tmp);
                ang = tmp.angle();
                profile.angle(e) = ang(1);
            end
        end
        % save flow
        function [fileList] = saveFlow(obj,outPort)
            data = reshape(obj.position_field_midline.d,obj.position_field_midline.oS);
            data = permute(data,[2 3 1]);
            data = reshape(data,[size(data,1) size(data,2)*size(data,3)]);
            % push the xy cordinates
            fileList = outPort.push('xyCoord.csv',data);
            % push the root tip
            fileList = outPort.push('rootTipObject.mat',obj);
            %
        end
    end
    
    
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save profile
        function [fileList] = save_flowProfile(profile,outPort)
           switch profile.type
                case 'steadyState'
                    fileList = rootTipSequence.save_steadyStateKinematicsProfile(profile,outPort);
                case 'nonSteadyState'
                    fileList = rootTipSequence.save_nonSteadyStateKinematicsProfile(profile,outPort);
           end
        end
        % plot flow profile
        function [h] = plot_flowProfile(profile)
           switch profile.type
                case 'steadyState'
                    h = rootTipSequence.plot_steadyStateKinematicsProfile(profile);
                case 'nonSteadyState'
                    h = rootTipSequence.plot_nonSteadyStateKinematicsProfile(profile);
           end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot steady state
        function [h] = plot_steadyStateKinematicsProfile(profile)
            h = figure;
            plot(profile.raw_domain,profile.raw_velocity,'.');
            hold on;
            [a g1 g2] = plotyy(profile.lDomain,profile.velocityProfile,profile.lDomain,profile.strainProfile);
            set(g1,'Color','r');
            
            growthRate = mean(profile.velocityProfile(end-10:end));
            maxY = 2*growthRate;
            
            maxStrain = max(profile.strainProfile);
            
            title({['Red-Percent Zone@' num2str(profile.PERZONE_WIDTH) '---' 'Cyan-Evan Zone@' num2str(profile.EVANSZONE_WIDTH)],['Yellow-Absolute Zone@' num2str(profile.ABSZONE_WIDTH)]});
            
            y1 = axis(a(1));
            y1(1) = 0;
            y1(2) = profile.lDomain(end);
            y1(3) = -1;
            y1(4) = maxY;            
            axis(a(1),y1);
            
            y1 = axis(a(2));
            y1(1) = 0;
            y1(2) = profile.lDomain(end);
            y1(3) = 0;
            y1(4) = maxStrain;    
            axis(a(2),y1);
            
            z(1) = axes;            
            plot(z(1),profile.lDomain,maxStrain*profile.PERZONE,'r');
            set(z(1),'Color','none');
            set(z(1),'Visible','off');
            y1 = axis(z(1));
            y1(1) = 0;
            y1(2) = profile.lDomain(end); 
            y1(3) = 0;
            y1(4) = maxStrain; 
            axis(z(1),y1);
            
           
            z(2) = axes;
            plot(z(2),profile.lDomain,maxStrain*profile.EVANSZONE,'c');
            set(z(2),'Color','none');
            set(z(2),'Visible','off');
            y1 = axis(z(2));
            y1(1) = 0;
            y1(2) = profile.lDomain(end);            
            y1(3) = 0;
            y1(4) = maxStrain; 
            axis(z(2),y1);
            
            z(3) = axes;
            plot(z(3),profile.lDomain,maxStrain*profile.ABSZONE,'y');
            set(z(3),'Color','none');
            set(z(3),'Visible','off');
            y1 = axis(z(3));
            y1(1) = 0;
            y1(2) = profile.lDomain(end);            
            y1(3) = 0;
            y1(4) = maxStrain; 
            axis(z(3),y1);
            
            
            xlabel('Along Root (px)');
            
            
            
        end
        % save profile
        function [fileList] = save_steadyStateKinematicsProfile(profile,outPort)
            h = rootTipSequence.plot_steadyStateKinematicsProfile(profile);
            fileList{1} = outPort.push('REGR.tif',h);
            fileList{2} = outPort.push('rawData.csv',[profile.raw_domain profile.raw_velocity]);
            fileList{3} = outPort.push('fitData.csv',[profile.lDomain' profile.velocityProfile' profile.strainProfile']);            
            % peak analysis
            fileList{4} = outPort.push('REGR_peak_location.csv',profile.maxLocation);
            fileList{5} = outPort.push('REGR_peak_value.csv',profile.maxREGR);
            % evans zone
            fileList{6} = outPort.push('EVANS_zone.csv',profile.EVANSZONE);
            fileList{7} = outPort.push('EVANS_zone_width.csv',profile.EVANSZONE_WIDTH);
            % percent zone
            fileList{8} = outPort.push('PERCENT_zone.csv',profile.PERZONE);
            fileList{9} = outPort.push('PERCENT_zone_width.csv',profile.PERZONE_WIDTH);
            % abs zone
            fileList{8} = outPort.push('ABS_zone.csv',profile.PERZONE);
            fileList{9} = outPort.push('ABS_zone_width.csv',profile.ABSZONE_WIDTH);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot steady state
        function [h] = plot_nonSteadyStateKinematicsProfile(profile)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make velocity figure
            h(1) = figure;
            mesh(profile.spaceDomain,profile.timeDomain,abs(profile.velocityProfile));
            view([0 90]);
            xlabel('Along Midline (px)');
            ylabel('Time (frames)');
            c1 = colorbar;
            caxis([0 6]);
            ylabel(c1,'Velocity (px/fr)');
            ax = axis;
            ax(4) = size(profile.spaceDomain,1);
            ax(3) = 1;
            axis(ax);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % REGR figure with nothing on it
            h(2) = figure;
            mesh(profile.spaceDomain,profile.timeDomain,abs(profile.strainProfile));
            view([0 90]);
            xlabel('Along Midline (px)');
            ylabel('Time (frames)');
            c2 = colorbar;
            caxis([0 7*10^-3]);
            ylabel(c2,'Strain (%/fr)');
            ax = axis;
            ax(4) = size(profile.spaceDomain,1);
            ax(3) = 1;
            axis(ax);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % REGR figure
            h(3) = figure;
            mesh(profile.spaceDomain,profile.timeDomain,abs(profile.strainProfile));
            view([0 90]);
            xlabel('Along Midline (px)');
            ylabel('Time (frames)');
            c2 = colorbar;
            caxis([0 7*10^-3]);
            ylabel(c2,'Strain (%/fr)');
            ax = axis;
            ax(4) = size(profile.spaceDomain,1);
            ax(3) = 1;
            axis(ax);
            
            % plot along the peak curve
            [val key] = max(profile.strainProfile,[],2);            
            tm = 1:numel(key);
            hold on
            plot3(key,tm,val,'k');
            
            Zone = [profile.profileSequence.EVANSZONE];
            hold on;
            % evans zone
            B = bwboundaries(Zone);
            for e = 1:numel(B)
                curve(e).Z = 2*Zone(B{e}(:,1),B{e}(:,2));
                curve(e).X = B{e}(:,1);
                curve(e).Y = B{e}(:,2);
                plot3(curve(e).X,curve(e).Y,curve(e).Z,'k','LineWidth',2,'Color','w','LineStyle','-');
            end
            
            Zone = [profile.profileSequence.PERZONE];
            hold on;
            % percent
            B = bwboundaries(Zone);
            for e = 1:numel(B)
                curve(e).Z = 2*Zone(B{e}(:,1),B{e}(:,2));
                curve(e).X = B{e}(:,1);
                curve(e).Y = B{e}(:,2);
                plot3(curve(e).X,curve(e).Y,curve(e).Z,'k','LineWidth',2,'Color','k','LineStyle','-');
            end
            
            Zone = [profile.profileSequence.ABSZONE];
            hold on;
            % abs percent
            B = bwboundaries(Zone);
            for e = 1:numel(B)
                curve(e).Z = 2*Zone(B{e}(:,1),B{e}(:,2));
                curve(e).X = B{e}(:,1);
                curve(e).Y = B{e}(:,2);
                plot3(curve(e).X,curve(e).Y,curve(e).Z,'k','LineWidth',2,'Color','y','LineStyle','-');
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % REGR figure with EVANS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % REGR figure
            h(4) = figure;
            mesh(profile.spaceDomain,profile.timeDomain,abs(profile.strainProfile));
            view([0 90]);
            xlabel('Along Midline (px)');
            ylabel('Time (frames)');
            c2 = colorbar;
            caxis([0 7*10^-3]);
            ylabel(c2,'Strain (%/fr)');
            ax = axis;
            ax(4) = size(profile.spaceDomain,1);
            ax(3) = 1;
            axis(ax);
            % plot along the peak curve
            [val key] = max(profile.strainProfile,[],2);            
            tm = 1:numel(key);
            hold on
            plot3(key,tm,val,'k');
            Zone = [profile.profileSequence.EVANSZONE];
            hold on;
            % evans zone
            B = bwboundaries(Zone);
            for e = 1:numel(B)
                curve(e).Z = 2*Zone(B{e}(:,1),B{e}(:,2));
                curve(e).X = B{e}(:,1);
                curve(e).Y = B{e}(:,2);
                plot3(curve(e).X,curve(e).Y,curve(e).Z,'k','LineWidth',2,'Color','w','LineStyle','-');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % REGR figure with EVANS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % REGR figure
            h(5) = figure;
            mesh(profile.spaceDomain,profile.timeDomain,abs(profile.strainProfile));
            view([0 90]);
            xlabel('Along Midline (px)');
            ylabel('Time (frames)');
            c2 = colorbar;
            caxis([0 7*10^-3]);
            ylabel(c2,'Strain (%/fr)');
            ax = axis;
            ax(4) = size(profile.spaceDomain,1);
            ax(3) = 1;
            axis(ax);
            % plot along the peak curve
            [val key] = max(profile.strainProfile,[],2);            
            tm = 1:numel(key);
            hold on
            plot3(key,tm,val,'k');
            Zone = [profile.profileSequence.PERZONE];
            hold on;
            % evans zone
            B = bwboundaries(Zone);
            for e = 1:numel(B)
                curve(e).Z = 2*Zone(B{e}(:,1),B{e}(:,2));
                curve(e).X = B{e}(:,1);
                curve(e).Y = B{e}(:,2);
                plot3(curve(e).X,curve(e).Y,curve(e).Z,'k','LineWidth',2,'Color','w','LineStyle','-');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % REGR figure with ABS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % REGR figure
            h(6) = figure;
            mesh(profile.spaceDomain,profile.timeDomain,abs(profile.strainProfile));
            view([0 90]);
            xlabel('Along Midline (px)');
            ylabel('Time (frames)');
            c2 = colorbar;
            caxis([0 7*10^-3]);
            ylabel(c2,'Strain (%/fr)');
            ax = axis;
            ax(4) = size(profile.spaceDomain,1);
            ax(3) = 1;
            axis(ax);
            % plot along the peak curve
            [val key] = max(profile.strainProfile,[],2);            
            tm = 1:numel(key);
            hold on
            plot3(key,tm,val,'k');
            Zone = [profile.profileSequence.ABSZONE];
            hold on;
            % evans zone
            B = bwboundaries(Zone);
            for e = 1:numel(B)
                curve(e).Z = 2*Zone(B{e}(:,1),B{e}(:,2));
                curve(e).X = B{e}(:,1);
                curve(e).Y = B{e}(:,2);
                plot3(curve(e).X,curve(e).Y,curve(e).Z,'k','LineWidth',2,'Color','w','LineStyle','-');
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % make velocity figure
            figure;
            growthRate = mean(profile.velocityProfile(:,end-10:end),2);
            hold on;
            plotyy(1:numel(profile.profileSequence),[profile.profileSequence.maxREGR],1:numel(profile.profileSequence),growthRate)
            
            
            
        end
        % save profile
        function [fileList] = save_nonSteadyStateKinematicsProfile(profile,outPort)
            h = rootTipSequence.plot_nonSteadyStateKinematicsProfile(profile);
            fileList = {};
            fileList{end+1} = outPort.push('velocity.tif',h(1));
            fileList{end+1} = outPort.push('REGR.tif',h(2));
            fileList{end+1} = outPort.push('REGR_labels.tif',h(3));
            fileList{end+1} = outPort.push('REGR_EVANS.tif',h(4));
            fileList{end+1} = outPort.push('REGR_PERCENT.tif',h(5));
            fileList{end+1} = outPort.push('REGR_ABS.tif',h(6));
            fileList{end+1} = outPort.push('strainProfile.csv',profile.strainProfile);
            fileList{end+1} = outPort.push('velocityProfile.csv',profile.velocityProfile);            
            % peak analysis
            % location
            fileList{end+1} = outPort.push('REGR_peak_location.csv',[profile.profileSequence.maxLocation]);
            fileList{end+1} = outPort.push('mean_REGR_peak_location.csv',mean([profile.profileSequence.maxLocation]));
            fileList{end+1} = outPort.push('std_REGR_peak_location.csv',std([profile.profileSequence.maxLocation]));
            % REGR value
            fileList{end+1} = outPort.push('REGR_peak_value.csv',[profile.profileSequence.maxREGR]);
            fileList{end+1} = outPort.push('mean_REGR_peak_value.csv',mean([profile.profileSequence.maxREGR]));
            fileList{end+1} = outPort.push('std_REGR_peak_value.csv',std([profile.profileSequence.maxREGR]));
            % evans zone
            fileList{end+1} = outPort.push('EVANS_zone.csv',[profile.profileSequence.EVANSZONE]);
            fileList{end+1} = outPort.push('EVANS_zone_width.csv',[profile.profileSequence.EVANSZONE_WIDTH]);
            fileList{end+1} = outPort.push('mean_EVANS_zone_width.csv',mean([profile.profileSequence.EVANSZONE_WIDTH]));
            fileList{end+1} = outPort.push('std_EVANS_zone_width.csv',std([profile.profileSequence.EVANSZONE_WIDTH]));
            % percent zone
            fileList{end+1} = outPort.push('PERCENT_zone.csv',[profile.profileSequence.PERZONE]);
            fileList{end+1} = outPort.push('PERCENT_zone_width.csv',[profile.profileSequence.PERZONE_WIDTH]);
            fileList{end+1} = outPort.push('mean_PERCENT_zone_width.csv',mean([profile.profileSequence.PERZONE_WIDTH]));
            fileList{end+1} = outPort.push('std_PERCENT_zone_width.csv',std([profile.profileSequence.PERZONE_WIDTH]));
            % abs zone
            fileList{end+1} = outPort.push('ABS_zone.csv',[profile.profileSequence.ABSZONE]);
            fileList{end+1} = outPort.push('ABS_zone_width.csv',[profile.profileSequence.ABSZONE_WIDTH]);
            fileList{end+1} = outPort.push('mean_ABS_zone_width.csv',mean([profile.profileSequence.ABSZONE_WIDTH]));
            fileList{end+1} = outPort.push('std_ABS_zone_width.csv',std([profile.profileSequence.ABSZONE_WIDTH]));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % fit strain
        function [modelFunction] = fitStrain(X,Y,modelType)
            modelFunction = modelType();
            modelFunction.fit(X,Y);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        % eval the model    
        function [Y,Yp] = evalStrainModel(model,Domain)
            %Y = veloc_SPEC(modelParameters,Domain);            
            Y = model.fnval(Domain);
            Yp = gradient(Y).*gradient(Domain).^-1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % generate metrics for flow profile
        function [profile] = generateStrainMetrics(profile,para)
            RELcutoff = para.RELcutoff;
            PERcutoff = para.PERcutoff;
            ABScutoff = para.ABScutoff;
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % analysis of width
            % get sig
            tmp = profile.strainProfile;                    
            % evans
            M = max(tmp);
            cRELcutoff = M - RELcutoff*M;
            profile.EVANSZONE = (tmp > cRELcutoff)';
            profile.EVANSZONE_WIDTH = sum(gradient(profile.lDomain).*profile.EVANSZONE');
            
            % percent zone
            [sdata sidx] = sort(tmp,'descend');
            sdata = cumsum(sdata);
            sdata = sdata/sdata(end);
            fidx = find(sdata < PERcutoff);
            z = zeros(size(sdata));
            z(sidx(fidx)) = 1;
            profile.PERZONE =  z';
            profile.PERZONE_WIDTH = sum(gradient(profile.lDomain).*profile.PERZONE');
            % abs zone
            profile.ABSZONE = (tmp > ABScutoff)';
            profile.ABSZONE_WIDTH = sum(gradient(profile.lDomain).*profile.ABSZONE');
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % analysis of peak and location of peak
            [profile.maxREGR profile.maxLocation] = max(tmp);
            profile.maxLocation = profile.lDomain(profile.maxLocation);
        end
    end
    
end