

%Generate simulation results
%Density_estimation_lastlast;

model_scale = 1;

layer_stl = 'squid.stl';
[Vertices_b,faces,normals,name] = stlRead(layer_stl);

layer_object.vertices = Vertices_b*model_scale;
layer_object.faces = faces;
layer_color = [1,0,0]; %Red squid

patch((layer_object),'FaceColor',       layer_color, ...
     'EdgeColor',       'none',        ...
     'FaceLighting',    'gouraud',     ...
     'AmbientStrength', 0.15);
set(gca,'Zdir','reverse','Ydir','reverse')
camlight('headlight');
material('shiny');
axis equal;
view(45, 35);
grid on;
