function [s] = PlotNMS3D(filename,PlotTitle)
  % PLOTNMS3D Display a NanoFocus NMS Dataset as 3-dimensional surface.
  %
  % Syntax:
  %   PlotNMS3D(filename,title)
  %
  % Parameters:
  %   filename - File specification for the NMS file to read and display
  %   title    - Plot title to display
  %
  % Example:
  %   PlotNMS3D('1-euro-star.nms','Star on a 1 Euro coin')

  % Copyright 2009, Dr. Georg Wiora, NanoFocus AG

  % Read in NMS-File
  s = openNMS(filename);
  % This returns a structure with the following fields:
  %   s.z       - M x N double matrix containing the depth values in micro meter
  %   s.u       - 1 x N vector containing the u-coordinates in micro meter for each row in z
  %   s.v       - M x 1 vector containing the v-coordinates in micro meter for each columnb in z
  %   s.r       - M x N matrix of UINT8 containing the reflection data. Values >0 denote valid pixels
  %   s.offsetu - double scalar containing the offset in micro meter for the u-coordinates
  %   s.offsetv - double scalar containing the offset in micro meter for the v-coordinates
  %   s.offsetz - double scalar containing the offset in micro meter for the z-coordinates
  %   s.comment - String array containing the files comment field
  %
  % Other possible call syntax is:
  % [u, v, z] = openNMS(filename)
  % [u, v, z, reflection] = openNMS(filename)

  % Display file comment if any
  fprintf(1,'File Comment: "%s"\n',s.comment);

  % Create new figure
  figure1 = figure('name',PlotTitle);

  % Create axes
  axes1 = axes('Parent',figure1,'FontWeight','bold','FontSize',14);
  view([-20 46]);
  hold('all');

  % Create surface
  surface('Parent',axes1,'ZData',s.z,'YData',s.v,'XData',s.u,...
    'LineStyle','none','CData',s.z);
  
  % Option 1: Make all axis equal scale to give a realistic representation
  % axis equal
  
  % Option 2: Make xy equal and overscale z by factor of 10 to show small height differences
  set(axes1,'DataAspectRatio',[1,1,1/10]);

  % Freeze axis ratio to allow nice 3d rotation
  axis vis3d

  % Create xlabel
  xlabel('X [�m]','FontWeight','bold','FontSize',14);

  % Create ylabel
  ylabel('Y [�m]','FontWeight','bold','FontSize',14);

  % Create zlabel
  zlabel('Z [�m]','FontWeight','bold','FontSize',14);

  % Create title
  title(['Profile ',PlotTitle],'Interpreter','none','FontWeight','bold',...
    'FontSize',16);

  % Create 2 lights to have a shiny surface. This works best with very smooth
  % surfaces. On rough surfaces it does not look very nice.
  % A spot light to create nice highlights
  light('Parent',axes1,'Style','local','Position',[-100, -400, 100]);
  % A dimmed spot light to provide a minimum brightness everywhere
  light('Parent',axes1,'Style','local','Position',[200, 200, 1000],'color',[0.5,0.5,0.5]);

  % Create colorbar
  colorbar('peer',axes1,'FontWeight','bold','FontSize',14);
  
  
  % Show reflection data
  figure('name',['Reflectivity ',PlotTitle]);
  imagesc(s.u,s.v,s.r);
  % Correct axis orientation
  axis xy
  colormap gray
  colorbar
  xlabel 'X [�m]'
  ylabel 'Y [�m]'
end
