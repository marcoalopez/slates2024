%%
clear variables

%% Set crystal symmetry and reference frame
CS = {... 
  'notIndexed',...
  crystalSymmetry('-3m1', [4.9 4.9 5.5], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Quartz', 'color', '#d9d9d9'),...
  crystalSymmetry('12/m1', [5.2 8.9 20], [90,95.8,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Muscovite', 'color', '#fb8072'),...
  crystalSymmetry('4/mmm', [4.6 4.6 3], 'mineral', 'Rutile', 'color', '#ffed6f'),...
  crystalSymmetry('-1', [8.1 13 7.2], [94.19,116.61,87.68]*degree, 'X||a*', 'Z||c', 'mineral', 'Low albite', 'color', '#80b1d3'),...
  crystalSymmetry('12/m1', [5.4 9.3 14], [90,96.282,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Chlorite', 'color', '#b3de69'),...
  crystalSymmetry('6/m', [9.4 9.4 6.9], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Hydroxylapatite', 'color', '#fdb462'),...
  crystalSymmetry('-3m1', [5 5 17], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Calcite', 'color', '#ccebc5'),...
  crystalSymmetry('m-3', [5.4 5.4 5.4], 'mineral', 'Pyrite', 'color', [0.58 0 0.83]),...
  crystalSymmetry('4/mmm', [6.6 6.6 6], 'mineral', 'Zircon', 'color', [0.5 0.5 0.5]),...
  crystalSymmetry('mmm', [5.8 3.5 6], 'mineral', 'Pyrrhotite high', 'color', [0.85 0.75 0.85]),...
  crystalSymmetry('-3', [5.1 5.1 14], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ilmenite', 'color', '#bc80bd')};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

%% Path to files
pname = 'C:\Users\Marco\_Documentos_\project_SLATES\raw_data\BEI\EBSD';
fname = [pname '\BEI_Map.ctf'];

%% Import Data
ebsd = EBSD.load(fname, CS, 'interface', 'ctf', 'convertEuler2SpatialReferenceFrame');

% rotate 90 degrees clockwise along the x vector to convert XZ to XY
ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, 90*degree), 'keepXY')

%% CLEAUP ROUTINE AND GRAIN RECONSTRUCTION
ebsd = ebsd('indexed');

% Remove wild spikes
[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle', 10*degree);
ebsd(grains(grains.grainSize < 2)) = [];

% Remove grain below 4 pixels
[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle', 10*degree);
ebsd(grains(grains.grainSize < 4)) = [];

% Reconstruct grains
[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle', 10*degree, 'boundary', 'tight');
grains


%% plot phase map with reconstructed grains
grains = smooth(grains, 2);

figure(1)
plot(ebsd)
hold on
plot(grains.boundary)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MUSCOVITE
% Define the Miller indices to plot in the pole figures
h = Miller({0,0,1}, {0,1,0}, ebsd('Muscovite').CS);

% Estimate ODF
psi = calcKernel(grains('Muscovite').meanOrientation)
odf_ms = calcDensity(grains('Muscovite').meanOrientation,...
                    'weigths', grains('Muscovite').area,...
                    'kernel', psi,...
                    'resolution', 1)

%% Pole figure
figure(2)
plotPDF(odf_ms, h, 'antipodal')
hold on
plotPDF(odf_ms, h, 'contour', [1,5,10,20,30], 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')


colormap(brewermap(256, 'Reds'));
CLim(gcm,'equal');
c = mtexColorbar('FontSize', 20);
c.Label.String = 'multiples of a uniform distribution';
c.Label.FontSize = 22;

%% IPF from foliation plane
figure(3)
plotIPDF(odf_ms, vector3d.Z)
hold on
plotIPDF(odf_ms, vector3d.Z, 'contour', [1,10,20,30], 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')
% hold on
% plot(h, 'upper', 'labeled', 'backgroundColor', 'w')

colormap(brewermap(256, 'Reds'));
c = mtexColorbar('FontSize', 20);
c.Label.String = 'multiples of a uniform distribution';
c.Label.FontSize = 22;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHLORITE
% Define the Miller indices to plot in the pole figures
h = Miller({0,0,1}, {0,1,0}, ebsd('Chlorite').CS);

% Estimate ODF
psi = calcKernel(grains('Chlorite').meanOrientation)
odf_chl = calcDensity(grains('Chlorite').meanOrientation,...
                     'weigths', grains('Chlorite').area,...
                     'kernel', psi,...
                     'resolution', 1)

%% Pole figure
figure(4)
plotPDF(odf_chl, h, 'antipodal')
hold on
plotPDF(odf_chl, h, 'contour', [1,5,20,40,60], 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')


colormap(brewermap(256, 'Greens'));
CLim(gcm,'equal');
c = mtexColorbar('FontSize', 20);
c.Label.String = 'multiples of a uniform distribution';
c.Label.FontSize = 22;

%% IPF from foliation plane
figure(5)
plotIPDF(odf_chl, vector3d.Z)
hold on
plotIPDF(odf_chl, vector3d.Z, 'contour', [1,5,20,40,60], 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')
% hold on
% plot(h, 'upper', 'labeled', 'backgroundColor', 'w')

colormap(brewermap(256, 'Greens'));
c = mtexColorbar('FontSize', 20);
c.Label.String = 'multiples of a uniform distribution';
c.Label.FontSize = 22;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUARTZ
% Define the Miller indices to plot in the pole figures
h = Miller({0,0,0,1}, {1,1,-2,0}, ebsd('Quartz').CS);

% Estimate ODF
psi = calcKernel(grains('Quartz').meanOrientation)
odf_Qtz = calcDensity(grains('Quartz').meanOrientation,...
                     'weigths', grains('Quartz').area,...
                     'kernel', psi,...
                     'resolution', 1)

%%
figure(6)
plotPDF(odf_Qtz, h, 'antipodal')
hold on
plotPDF(odf_Qtz, h, 'contour', 1:0.2:2, 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')


%colormap(flipud(colormap('bone')));
colormap(brewermap(256, 'GnBu'));
CLim(gcm,'equal');
c = mtexColorbar('FontSize', 20);
c.Label.String = 'multiples of a uniform distribution';
c.Label.FontSize = 22;

%%
figure(7)
plotIPDF(odf_Qtz, vector3d.Z)
hold on
plotIPDF(odf_Qtz, vector3d.Z, 'contour', 1:0.3:2, 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')

%colormap(flipud(colormap('bone')))
colormap(brewermap(256, 'GnBu'));
c = mtexColorbar('FontSize', 20);
c.Label.String = 'multiples of a uniform distribution';
c.Label.FontSize = 22;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLAGIOCLASE
% Define the Miller indices to plot in the pole figures
h = Miller({1,0,0}, {0,1,0}, {0,0,1}, ebsd('Low albite').CS);

% Estimate ODF
psi = calcKernel(grains('Low albite').meanOrientation)
odf_Pl = calcDensity(grains('Low albite').meanOrientation,...
                    'weigths', grains('Low albite').area,...
                    'kernel', psi,...
                     'resolution', 1)

%% Pole figure
figure(8)
plotPDF(odf_Pl, h, 'antipodal')
hold on
plotPDF(odf_Pl, h, 'contour', 1:0.5:3, 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')


colormap(brewermap(256, 'Blues'));
CLim(gcm,'equal');
c = mtexColorbar('FontSize', 20);
c.Label.String = 'multiples of a uniform distribution';
c.Label.FontSize = 22;

%% IPF from foliation plane
figure(9)
plotIPDF(odf_Pl, vector3d.Z)
hold on
plotIPDF(odf_Pl, vector3d.Z, 'contour', 1:0.5:3, 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')
hold on
plot(h, 'upper', 'labeled', 'backgroundColor', 'w','MarkerFaceColor','Gold','MarkerEdgeColor','black')

colormap(brewermap(256, 'Blues'));
c = mtexColorbar('FontSize', 20);
c.Label.String = 'multiples of a uniform distribution';
c.Label.FontSize = 22;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate CPO/texture indexes (J index) from main mineral phases
j_ms = norm(odf_ms)^2;
j_chl = norm(odf_chl)^2;
j_plag = norm(odf_Pl)^2;
j_qtz = norm(odf_Qtz)^2;

%%
disp ' ';
disp 'J (texture) index values'
fprintf('Muscovite = %.3f \n', j_ms);
fprintf('Chlorite = %.3f \n', j_chl);
fprintf('Quartz = %.3f \n', j_qtz);
fprintf('Plagioclase = %.3f \n', j_plag);
disp ' ';

%%
e_ms = entropy(odf_ms);
e_chl = entropy(odf_chl);

disp ' ';
disp 'ODF entropy'
fprintf('Muscovite = %.3f \n', e_ms);
fprintf('Chlorite = %.3f \n', e_chl);
disp ' ';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate the number of the main mineral phases (needed to asses the quality of the CPO and J index)
% comment: TODO
muscovite = calcGrains(ebsd('Muscovite'), 'angle', 10*degree);
chlorite= calcGrains(ebsd('Chlorite'), 'angle', 10*degree);
plagiclase = calcGrains(ebsd('Low albite'), 'angle', 10*degree);
quartz = calcGrains(ebsd('Quartz'), 'angle', 10*degree);

%%
disp ' ';
disp 'NUMBER OF GRAINS'
fprintf('Muscovite = %.0f \n', length(muscovite.id));
fprintf('Chlorite = %.0f \n', length(chlorite.id));
fprintf('Quartz = %.0f \n', length(quartz.id));
fprintf('Plagioclase = %.0f \n', length(plagiclase.id));
disp ' '

%% IPF map (foliation plane)
% % Set IPF colors respect to foliation plane (yvector)
% oM = ipfHSVKey(ebsd('Muscovite'));
% oM.inversePoleFigureDirection = zvector;
% IPF_z = oM.orientation2color(ebsd('Muscovite').orientations);
% 
% % plot IPF map
% figure(4)
% plot(ebsd('Muscovite'), IPF_z);
% hold on
% plot(grains('Muscovite').boundary, 'linewidth', 0.5, 'linecolor', [0.5 0.5 0.5], 'FaceAlpha', 0.5)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST ODF harmonic vs ODF kde vs  raw

% % CHLORITE
% % Define the Miller indices to plot in the pole figures
% h = Miller({0,0,1}, ebsd('Chlorite').CS);
% 
% % Estimate ODF
% psi = calcKernel(grains('Chlorite').meanOrientation)
% odf_chl2 = calcKernelODF(grains('Chlorite').meanOrientation,...
%                        'weigths', grains('Chlorite').area,...
%                        'kernel', psi)
% 
% %% plot
% mtexFig = newMtexFigure('figSize', 'huge', 'layout', [1 3]);
% 
% mtexTitle('ODF harmonic method', 'doNotDraw')
% plotPDF(odf_chl, h, 'antipodal')
% hold on
% plotPDF(odf_chl, h, 'contour', [10,20,30,40,50,60], 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')
% colormap(brewermap(256, 'Greens'));
% c = mtexColorbar('FontSize', 20);
% c.Label.String = 'multiples of a uniform distribution';
% c.Label.FontSize = 22;
% 
% nextAxis(1,2)
% mtexTitle('ODF KDE method', 'doNotDraw')
% plotPDF(odf_chl2, h, 'antipodal')
% hold on
% plotPDF(odf_chl2, h, 'contour', [10,20,30,40,50,60], 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')
% colormap(brewermap(256, 'Greens'));
% c = mtexColorbar('FontSize', 20);
% c.Label.String = 'multiples of a uniform distribution';
% c.Label.FontSize = 22;
% 
% nextAxis(1,3)
% mtexTitle('Raw data', 'doNotDraw')
% plotPDF(grains('Chlorite').meanOrientation, h, 'weights', grains('Chlorite').area,...
%         'antipodal', 'contourf', [10,20,30,40,50,60], 'linecolor', 'black', 'linewidth', 2, 'ShowText', 'on')
% colormap(brewermap(256, 'Greens'));
% c = mtexColorbar('FontSize', 20);
% 
% %%
% norm(odf_chl2)^2
