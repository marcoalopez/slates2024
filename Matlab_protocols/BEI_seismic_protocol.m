%% Sample reference: BEI
% *Description*
% 
% *TODO...*micas represent the 47.50 % of the rock and Chl/Ms ration is 0.43
% 
% *Mineral content* (estimated using EDX chemical maps)
% 
% Muscovite = 33.2 %
% 
% Chlorite = 14.3 %
% 
% Quartz = 24.0 %
% 
% Plagioclase = 21.8 %
% 
% Others = 6 %
%% Import EBSD data

%%
clear variables

% Set crystal symmetry and reference frame
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

% Import Data
pname = 'C:\Users\Marco\_Documentos_\project_SLATES\raw_data\BEI\EBSD';
fname = [pname '\BEI_Map.ctf'];
ebsd = EBSD.load(fname, CS, 'interface', 'ctf', 'convertEuler2SpatialReferenceFrame');

% rotate 90 degrees clockwise along the x vector to convert XZ to XY
ebsd = rotate(ebsd, rotation.byAxisAngle(xvector, 90*degree), 'keepXY')
%% Cleanup routine and grain reconstruction

ebsd = ebsd('indexed');

% Remove wild spikes
[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle', 10*degree);
ebsd(grains(grains.grainSize < 2)) = [];

% Remove grains below 4 pixels
[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle', 10*degree);
ebsd(grains(grains.grainSize < 4)) = [];

% Reconstruct grains
[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, 'angle', 10*degree, 'boundary', 'tight');
grains
%% Estimate orientation distribution functions (ODFs) for the main mineral phases
% TODO

% MUSCOVITE
psi = calcKernel(grains('Muscovite').meanOrientation)
odf_ms = calcKernelODF(grains('Muscovite').meanOrientation,...
                       'weigths', grains('Muscovite').area,...
                       'kernel', psi)

% CHLORITE
psi = calcKernel(grains('Chlorite').meanOrientation)
odf_chl = calcKernelODF(grains('Chlorite').meanOrientation,...
                        'weigths', grains('Chlorite').area,...
                        'kernel', psi)

% QUARTZ
psi = calcKernel(grains('Quartz').meanOrientation)
odf_qtz = calcKernelODF(grains('Quartz').meanOrientation,...
                        'weigths', grains('Quartz').area,...
                        'kernel', psi)

% PLAGIOCLASE
psi = calcKernel(grains('Low albite').meanOrientation)
odf_alb = calcKernelODF(grains('Low albite').meanOrientation,...
                        'weigths', grains('Low albite').area,...
                        'kernel', psi)
%% Define the elastic stiffness tensors of major phases
% Alpha quartz
% We will use the elastic constants and crystal lattice paremeters provided 
% by
%% 
% * Wang, J., Mao, Z., Jiang, F., Duffy, T.S., 2015. Elasticity of single-crystal 
% quartz to 10 GPa. Phys Chem Minerals 42, 203–212. <https://doi.org/10.1007/s00269-014-0711-z 
% https://doi.org/10.1007/s00269-014-0711-z>
%% 
% at room pressure.

qtz_density = 2.648;  % (g/cm3)

% Reference frame for elastic constants at 1 atm
cs_quartz = crystalSymmetry('-3m1',...
                             [4.913, 4.913, 5.405],...
                             [90.0, 90.0, 120.0]*degree,...
                             'x||a', 'z||c', 'mineral', 'alpha_Quartz');

% independent elastic constants of alpha quartz (trigonal = 6) at RP-RT
C11 =  86.6;
C33 = 106.4;
C44 =  58.0;
C12 =  6.74;
C13 =  12.4;
C14 =  17.8;

% dependent terms
C66 = 0.5 * (C11 - C12);
C22 = C11; C55 = C44; C23 = C13; C24 = -C14; C56 = C14;


% Elastic stiffness tensor (in GPa) values as a matrix Cij (determined at
% 30 C)
Cij_aQtz2006 =... 
    [[C11   C12   C13   C14   0.0   0.0];...
    [ C12   C22   C23   C24   0.0   0.0];...
    [ C13   C13   C33   0.0   0.0   0.0];...
    [ C14   C24   0.0   C44   0.0   0.0];...
    [ 0.0   0.0   0.0   0.0   C55   C56];...
    [ 0.0   0.0   0.0   0.0   C56   C66]];

% Determine the elastic stiffness tensor (Cij) of alpha quartz Units are in 'GPa'
tensor_quartz = stiffnessTensor(Cij_aQtz2006,...
                                'unit', 'GPa',...
                                cs_quartz,...
                                'density', qtz_density);
% Low albite
% We will use the elastic constants and crystal lattice paremeters provided 
% in::
%% 
% * J.M. Brown, E.H. Abramson and R.J. Angel (2006) Triclinic elastic constants 
% for low albite. _Phys. Chem. Minerals_ 33: 256-265 <https://doi.org/10.1007/s00269-006-0074-1 
% https://doi.org/10.1007/s00269-006-0074-1> 
%% 
% In this, the Y tensor axis is aligned parallel to crystal b direction of the 
% real lattice (Y || b). The X tensor axis is aligned parallel to crystal a* direction 
% of the reciprocal lattice (X || a*). The a* direction is perpendicular the b 
% and c directions of the direct lattice of all crystals of any symmetry, including 
% triclinic.

alb_density = 2.623;  % source Almqvist et al. (2017), (g/cm3)

% Non-standard crystal reference frame for elastic constants of triclinic albite crystal
cs_alb2006 = crystalSymmetry('-1', [8.13662, 12.7857, 7.1582],...
                                    [94.253, 116.605, 87.756]*degree,...
                                    'X||a*', 'Z||c', 'mineral', 'Albite-2006');

% independent elastic constants of low albite (triclinic = 21) at RP-RT
C11 =  69.1;
C22 = 183.5;
C33 = 179.5;
C44 =  24.9;
C55 =  26.8;
C66 =  33.5;
C12 =  34.0;
C13 =  30.8;
C14 =   5.1;
C15 =  -2.4;
C16 =  -0.9;
C23 =   5.5;
C24 =  -3.9;
C25 =  -7.7;
C26 =  -5.8;
C34 =  -8.7;
C35 =   7.1;
C36 =  -9.8;
C45 =  -2.4;
C46 =  -7.2;
C56 =   0.5;

% Elastic stiffness tensor (in GPa) values as a matrix Mij
Cij_alb2006 =...
    [[C11   C12   C13   C14   C15   C16];...
    [ C12   C22   C23   C24   C25   C26];...
    [ C13   C23   C33   C34   C35   C36];...
    [ C14   C24   C34   C44   C45   C46];...
    [ C15   C25   C35   C45   C55   C56];...
    [ C16   C26   C36   C46   C56   C66]];

% Determine the elastic stiffness tensor (Cij) of alpha quartz Units are in 'GPa'
tensor_albite = stiffnessTensor(Cij_alb2006,...
                                'unit', 'GPa',... 
                                cs_alb2006,...
                                'density', alb_density);
% Muscovite
% We will use the elastic constants and crystal lattice paremeters provided 
% in:
%% 
% * Militzer, B., Wenk, H.-R., Stackhouse, S., Stixrude, L., 2011. First-principles 
% calculation of the elastic moduli of sheet silicates and their application to 
% shale anisotropy. American Mineralogist 96, 125–137. <https://doi.org/10.2138/am.2011.3558 
% https://doi.org/10.2138/am.2011.3558>
% * M.T. Vaughan and S. Guggenheim (1986) Elasticity of Muscovite and its relationship 
% to crystal structure. _J. Gephys. Res._ Vol.91 pp.4657-4664. <https://doi.org/10.1029/JB091iB05p04657 
% https://doi.org/10.1029/JB091iB05p04657>

ms_density = 2.830;  % source Almqvist et al. (2017), in g/cm3

% Reference frame for elastic constants
cs_Ms1986 = crystalSymmetry('12/m1', [5.1579, 8.9505, 20.0710],...
                                  [90.0, 95.75, 90.0]*degree,...
                                  'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Muscovite-1986');

% independent elastic constants of muscovite (monoclinic = 13) at RP-RT
C11 = 180.9;
C22 = 170.0;
C33 =  60.3;
C44 =  18.4;
C55 =  23.8;
C66 =  70.5;
C12 =  53.4;
C13 =  27.2;
C15 = -14.7;
C23 =  23.5;
C25 =   1.4;
C35 =  -1.0;
C46 =  -1.8;

% Elastic stiffness tensor (in GPa) values as a matrix Mij (IRE convention)
Cij_Ms2011 =...
    [[C11   C12   C13   0.0   C15   0.0];...
    [ C12   C22   C23   0.0   C25   0.0];...
    [ C13   C23   C33   0.0   C35   0.0];...
    [ 0.0   0.0   0.0   C44   0.0   C46];...
    [ C15   C25   C35   0.0   C55   0.0];...
    [ 0.0   0.0   0.0   C46   0.0   C66]];


% Determine the elastic stiffness tensor (Cij) of alpha quartz Units are in 'GPa'.
tensor_muscovite = stiffnessTensor(Cij_Ms2011,...
                                   'unit', 'GPa',...
                                   cs_Ms1986,...
                                   'density', ms_density);
% Clinochlore (Chlorite)
% We will use the elastic constants and crystal lattice paremeters provided 
% in:
%% 
% * Joswig, W., Feuss, H. and Mason, S.A. (1989) Neutron diffraction study of 
% a one-layer monoclinic chlorite. _Clays and Clay Miner_. 37:511-514. <https://doi.org/10.1346/CCMN.1989.0370602 
% https://doi.org/10.1346/CCMN.1989.0370602> 
% * Mookherjee, M., Mainprice, D., 2014. Unusually large shear wave anisotropy 
% for chlorite in subduction zone settings. Geophysical Research Letters 41, 1506–1513. 
% <https://doi.org/10.1002/2014GL059334 https://doi.org/10.1002/2014GL059334>

chl_density = 2.71;  % averaged from several chlinochlores in Katahara (1996), in g/cm3

% Reference frame for elastic constants
cs_chl1989 = crystalSymmetry('12/m1', [5.327, 9.227, 14.327],...
                            [90.00, 96.81, 90.00]*degree,...
                            'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Clinochlore IIb-2 1989');

% Table 3 in Mookherjee and Mainprice (2014) RP-RT
C11 = 197.8;
C22 = 202.3;
C33 = 135.1;
C44 =  24.5;
C55 =  24.4;
C66 =  70.3;
C12 =  60.7;
C13 =  21.1;
C15 =   3.3;
C23 =  34.1;
C25 =   0.2;
C35 =   0.4;
C46 =   0.1;

% Elastic stiffness tensor (in GPa) values as a matrix Mij
Cij_Chl2014 =...
    [[C11   C12   C13   0.0   C15   0.0];...
    [ C12   C22   C23   0.0   C25   0.0];...
    [ C13   C23   C33   0.0   C35   0.0];...
    [ 0.0   0.0   0.0   C44   0.0   C46];...
    [ C15   C25   C35   0.0   C55   0.0];...
    [ 0.0   0.0   0.0   C46   0.0   C66]];

% Determine the elastic stiffness tensor (Cij) of alpha quartz Units are in 'GPa'.
tensor_chlorite = stiffnessTensor(Cij_Chl2014,...
                                 'unit', 'GPa',...
                                 cs_chl1989,...
                                 'density', chl_density);
%% Estimate the rock stiffness/elastic tensor
% For this...TODO
% 
% The indexation of the different main mineral phases was notably different, 
% especially for the micas. Due to this, we have chosen to estimate the modal 
% proportions using EDX map data and image segmentation instead of using EBSD 
% data.
% 
% 

% Estimate the ODF-weighted elastic tensor for each mineral phase
qtz_weighted = calcTensor(odf_qtz, tensor_quartz, 'Hill');
alb_weighted = calcTensor(odf_alb, tensor_albite, 'Hill');
ms_weighted  = calcTensor(odf_ms, tensor_muscovite, 'Hill');
chl_weighted = calcTensor(odf_chl, tensor_chlorite, 'Hill');

% set the volume mineral content of main phases (recalculated to 100 %)
qtz_wt = 0.257;
alb_wt = 0.234;
ms_wt = 0.356;
chl_wt = 0.153;

% estimate the average elastic tensor of the rock
[C_voigt, C_reuss, C_hill] = mean([qtz_weighted, alb_weighted, ms_weighted, chl_weighted],...
                             'weights', [qtz_wt, alb_wt, ms_wt, chl_wt])
%%
% check density
rock_density = qtz_density * qtz_wt +...
               alb_density * alb_wt +...
               ms_density * ms_wt +...
               chl_density * chl_wt
%% Estimate the elastic wave velocities (km/s) and anisotropy of rocks
% Once we estimated the elastic stiffness tensor fo the rock/aggregate we need 
% to solve the Christoffel equation to estimate the wave velocities as a function 
% of the propagation direction and, thus, the anisotropy. For this we use the 
% |velocity| function that computes the elastic wave velocity (km/s). The inputs 
% are the elastic stiffness tensor (GPa) and the material density (g/cm3). As 
% an optional input, we can define a specific propagation direction or list of 
% propagation directions (@vector3d). By default, the funtion estimate the wave 
% velocities for all propagation directions. The function returns the following 
% estimates:
%% 
% * vp  - velocity of the p-wave (km/s)
% * vs1 - velocity of the s1-wave (km/s)
% * vs2 - velocity of the s2-wave ( km/s)
% * pp  - polarisation of the p-wave (particle movement, vibration direction)
% * ps1 - polarisation of the s1-wave (particle movement, vibration direction)
% * ps2 - polarisation of the s2-wave (particle movement, vibration direction)
%% 
% *Approaches*
% 
% *Voigt*: Determined using the stiffness tensor _Cij:_ provides the upper bound 
% for wave speeds
% 
% *Reuss*: Determined using the compliance tensor _Sij:_ provides the lower 
% bound for wave speeds
% 
% *Hill or VRH*: Average of the Voigt and Reuss bounds
% 
% 
% 
% *Anisotropy* (in percentage) is calculated as follows
% 
% 200 × (Vmax - Vmin) / (Vmax + Vmin)
% 
% while for shear wave splitting
% 
% 100 × Vs1i − Vs2i / mean(Vs1)

% Estimate elastic wave velocities (km/s) from the stiffness tensor
[vp, vs1, vs2, pp, ps1, ps2] = velocity(C_hill, 'harmonic');

% P-WAVES
vp_xyz = eval(vp, [xvector, yvector, zvector]);  % principal directions
[maxVp, maxVpPos] = max(vp);                     % max
[minVp, minVpPos] = min(vp);                     % min
AVp = 200 * (maxVp - minVp) ./ (maxVp + minVp);  % anisotropy (in percentage)

% S1-WAVES
vs1_xyz = eval(vs1, [xvector, yvector, zvector]);  % principal directions
[maxS1, maxS1pos] = max(vs1);                      % max
[minS1, minS1pos] = min(vs1);                      % min
AVs1 = 200  * (maxS1 - minS1) ./ (maxS1 + minS1);  % anisotropy (in percentage)

% S2-WAVES
vs2_xyz = eval(vs2, [xvector, yvector, zvector]);  % principal directions
[maxS2, maxS2pos] = max(vs2);                      % max
[minS2, minS2pos] = min(vs2);                      % min
AVs2 = 200  * (maxS2 - minS2) ./ (maxS2 + minS2);  % anisotropy (in percentage)

% Shear wave splitting: Velocity difference Vs1-Vs2 (km/s)
dVs = vs1 - vs2;
[maxdVs, maxdVsPos] = max(dVs);          % max
[mindVs, mindVsPos] = min(dVs);          % min
AVs = 200 * (vs1 - vs2) ./ (vs1 + vs2);  % S-wave anisotropy (in percentage)
[maxAVs, maxAVsPos] = max(AVs);          % max S-wave anisotropy (in percentage)
[minAVs, minAVsPos] = min(AVs);          % min S-wave anisotropy (in percentage)

% RATIOS
ratio_vpvs = vp ./ vs1;                                        % Vp/Vs1 ratio
[maxVpVs1, maxVpVs1Pos] = max(ratio_vpvs);                     % max
[minVpVs1, minVpVs1Pos] = min(ratio_vpvs);                     % min
AVpVs1 = 200 * (maxVpVs1 - minVpVs1) / (maxVpVs1 + minVpVs1);  % anisotropy (in percentage)

ratio_vpvs2 = vp ./ vs2;                                       % Vp/Vs2 ratio
[maxVpVs2, maxVpVs2Pos] = max(ratio_vpvs2);                    % max
[minVpVs2, minVpVs2Pos] = min(ratio_vpvs);                     % min
AVpVs2 = 200 * (maxVpVs2 - minVpVs2) / (maxVpVs2 + minVpVs2);  % anisotropy (in percentage)
%% Summarize results

disp ' ';
disp 'RESULTS:';
disp ' ';

disp 'P-WAVES'
fprintf('Maximum P-wave velocity = %.3f (km/s).\n', maxVp);
fprintf('Minimum P-wave velocity = %.3f (km/s).\n', minVp);
fprintf('P-wave anisotropy = %.1f (in percentage).\n', AVp);
disp '====================================================';

disp 'S-WAVES'
fprintf('Max. Vs1 = %.4f (km/s).\n', maxS1);
fprintf('Min. Vs1 = %.4f (km/s).\n', minS1);
fprintf('S1-wave anisotropy = %.1f (km/s).\n', AVs1);
disp '---';
fprintf('Max. Vs2 = %.4f (km/s).\n', maxS2);
fprintf('Min. Vs2 = %.4f (km/s).\n', minS2);
fprintf('S2-wave anisotropy = %.1f (km/s).\n', AVs2);
disp '====================================================';

disp 'Shear Wave Splitting'
fprintf('Max. velocity difference Vs1-Vs2 = %.4f (km/s).\n', maxdVs);
fprintf('Min. velocity difference Vs1-Vs2 = %.4f (km/s).\n', mindVs);
fprintf('SWS max. anisotropy = %.1f (in percentage).\n', maxAVs);
fprintf('SWS min. anisotropy = %.1f (in percentage).\n', minAVs);
disp '====================================================';

disp 'Vp/Vs1'
fprintf('Max. Vp/Vs1 = %.4f\n', maxVpVs1);
fprintf('Min. Vp/Vs1 = %.4f\n', minVpVs1);
fprintf('Vp/Vs1 Anisotropy = %.1f (in percentage).\n', AVpVs1);
disp '---';
disp 'Vp/Vs2'
fprintf('Max. Vp/Vs1 = %.4f\n', maxVpVs2);
fprintf('Min. Vp/Vs1 = %.4f\n', minVpVs2);
fprintf('Vp/Vs2 Anisotropy = %.1f (in percentage).\n', AVpVs2);
disp '====================================================';

disp ' ';
disp 'SEISMIC VELOCITIES PERPENDICULAR TO FOLIATION'
fprintf('P-wave velocity = %.3f (km/s)\n', vp_xyz(3));
fprintf('S1-wave velocity = %.3f (km/s)\n', vs1_xyz(3));
fprintf('S2-wave velocity = %.3f (km/s)\n', vs2_xyz(3));
fprintf('Velocity difference Vs1-Vs2 = %.4f (km/s)\n', vs1_xyz(3) - vs2_xyz(3));
fprintf('Vp/Vs1 = %.2f\n', vp_xyz(3) / vs1_xyz(3));
fprintf('Vp/Vs2 = %.2f\n', vp_xyz(3) / vs2_xyz(3));
disp ' ';
disp 'SEISMIC VELOCITIES PERPENDICULAR TO SURFACE MEASURED AND PARALLEL TO FOLIATION'
fprintf('P-wave velocity = %.3f (km/s)\n', vp_xyz(2));
fprintf('S1-wave velocity = %.3f (km/s)\n', vs1_xyz(2));
fprintf('S2-wave velocity = %.3f (km/s)\n', vs2_xyz(2));
fprintf('Velocity difference Vs1-Vs2 = %.4f (km/s)\n', vs1_xyz(2) - vs2_xyz(2));
fprintf('Vp/Vs1 = %.2f\n', vp_xyz(2) / vs1_xyz(2));
fprintf('Vp/Vs2 = %.2f\n', vp_xyz(2) / vs2_xyz(2))
disp ' ';
disp 'SEISMIC VELOCITIES PARALLEL TO SURFACE MEASURED AND THE FOLIATION'
fprintf('P-wave velocity = %.3f (km/s)\n', vp_xyz(1));
fprintf('S1-wave velocity = %.3f (km/s)\n', vs1_xyz(1));
fprintf('S2-wave velocity = %.3f (km/s)\n', vs2_xyz(1));
fprintf('Velocity difference Vs1-Vs2 = %.4f (km/s)\n', vs1_xyz(1) - vs2_xyz(1));
fprintf('Vp/Vs1 = %.2f\n', vp_xyz(1) / vs1_xyz(1));
fprintf('Vp/Vs2 = %.2f\n', vp_xyz(1) / vs2_xyz(1))
disp ' ';
%%
plotSeismicVelocities(C_hill)
colormap(flipud(brewermap(256, 'Spectral')))
% Estimate the ODF-weighted elastic tensor for each mineral phase
qtz_weighted = calcTensor(odf_qtz, tensor_quartz, 'geometric');
alb_weighted = calcTensor(odf_alb, tensor_albite, 'geometric');
ms_weighted  = calcTensor(odf_ms, tensor_muscovite, 'geometric');
chl_weighted = calcTensor(odf_chl, tensor_chlorite, 'geometric');

% set the volume mineral content of main phases (recalculated to 100 %)
qtz_wt = 0.257;
alb_wt = 0.234;
ms_wt = 0.356;
chl_wt = 0.153;

% estimate the average elastic tensor of the rock
C_geometric = mean([qtz_weighted, alb_weighted, ms_weighted, chl_weighted],...
                   'weights', [qtz_wt, alb_wt, ms_wt, chl_wt], 'geometric')
%%
% Estimate elastic wave velocities (km/s) from the stiffness tensor
[vp, vs1, vs2, pp, ps1, ps2] = velocity(C_geometric, 'harmonic');

% P-WAVES
vp_xyz = eval(vp, [xvector, yvector, zvector]);  % principal directions
[maxVp, maxVpPos] = max(vp);                     % max
[minVp, minVpPos] = min(vp);                     % min
AVp = 200 * (maxVp - minVp) ./ (maxVp + minVp);  % anisotropy (in percentage)

% S1-WAVES
vs1_xyz = eval(vs1, [xvector, yvector, zvector]);  % principal directions
[maxS1, maxS1pos] = max(vs1);                      % max
[minS1, minS1pos] = min(vs1);                      % min
AVs1 = 200  * (maxS1 - minS1) ./ (maxS1 + minS1);  % anisotropy (in percentage)

% S2-WAVES
vs2_xyz = eval(vs2, [xvector, yvector, zvector]);  % principal directions
[maxS2, maxS2pos] = max(vs2);                      % max
[minS2, minS2pos] = min(vs2);                      % min
AVs2 = 200  * (maxS2 - minS2) ./ (maxS2 + minS2);  % anisotropy (in percentage)

% Shear wave splitting: Velocity difference Vs1-Vs2 (km/s)
dVs = vs1 - vs2;
[maxdVs, maxdVsPos] = max(dVs);          % max
[mindVs, mindVsPos] = min(dVs);          % min
AVs = 200 * (vs1 - vs2) ./ (vs1 + vs2);  % S-wave anisotropy (in percentage)
[maxAVs, maxAVsPos] = max(AVs);          % max S-wave anisotropy (in percentage)
[minAVs, minAVsPos] = min(AVs);          % min S-wave anisotropy (in percentage)

% RATIOS
ratio_vpvs = vp ./ vs1;                                        % Vp/Vs1 ratio
[maxVpVs1, maxVpVs1Pos] = max(ratio_vpvs);                     % max
[minVpVs1, minVpVs1Pos] = min(ratio_vpvs);                     % min
AVpVs1 = 200 * (maxVpVs1 - minVpVs1) / (maxVpVs1 + minVpVs1);  % anisotropy (in percentage)

ratio_vpvs2 = vp ./ vs2;                                       % Vp/Vs2 ratio
[maxVpVs2, maxVpVs2Pos] = max(ratio_vpvs2);                    % max
[minVpVs2, minVpVs2Pos] = min(ratio_vpvs);                     % min
AVpVs2 = 200 * (maxVpVs2 - minVpVs2) / (maxVpVs2 + minVpVs2);  % anisotropy (in percentage)
%%
fprintf('Maximum P-wave velocity = %.3f (km/s).\n', maxVp);
%%
fprintf(['P-wave velocity (direction-z)= %.3f (km/s)\n',...
         'P-wave velocity (direction-y)= %.3f (km/s)\n',...
         'P-wave velocity (direction-x)= %.3f (km/s)\n'],...
         vp_xyz(3), vp_xyz(2), vp_xyz(1));
%%
fprintf(['S1-wave velocity (direction-z)= %.3f (km/s)\n',...
         'S1-wave velocity (direction-y)= %.3f (km/s)\n',...
         'S1-wave velocity (direction-x)= %.3f (km/s)\n'],...
         vs1_xyz(3), vs1_xyz(2), vs1_xyz(1));
%%
plotSeismicVelocities(C_geometric)
colormap(flipud(brewermap(256, 'Spectral')))
%%
% set a grid with a resolution of 1 degree
XY_grid = equispacedS2Grid('upper', 'resolution', 1*degree);

% estimate velocities using grid
[vp_hill, vs1_hill, vs2_hill, pp, ps1, ps2] = velocity(C_hill, XY_grid);
[vp_voigt, vs1_voigt, vs2_voigt, pp, ps1, ps2] = velocity(C_voigt, XY_grid);
[vp_reuss, vs1_reuss, vs2_reuss, pp, ps1, ps2] = velocity(C_reuss, XY_grid);
[vp_geo, vs1_geo, vs2_geo, pp, ps1, ps2] = velocity(C_geometric, XY_grid);

% create a table with Vp, and directions in polar coordinates
dataset = table(vp_voigt.',...
                vp_reuss.',...
                vp_hill.',...
                vp_geo.',...
                vs1_voigt.',...
                vs1_reuss.',...
                vs1_hill.',...
                vs1_geo.',...
                vs2_voigt.',...
                vs2_reuss.',...
                vs2_hill.',...
                vs2_geo.',... 
                XY_grid.theta.',...
                XY_grid.rho.',...
                'VariableNames',{'Vp_voigt', 'Vp_reuss', 'Vp_hill', 'Vp_geo',...
                                 'Vs1_voigt', 'Vs1_reuss', 'Vs1_hill', 'Vs1_geo',...
                                 'Vs2_voigt', 'Vs2_reuss', 'Vs2_hill', 'Vs2_geo',...
                                 'theta_rad','rho_rad'});

writetable(dataset, 'BEI_calc_velocities.csv', 'Delimiter', ';')