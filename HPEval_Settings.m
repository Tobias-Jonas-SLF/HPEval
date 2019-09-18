function hpeset            = HPEval_Settings

  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  % GENERAL DESCRIPTION                                                  %
  %  %   HPEval_Settings defines setting for HPEval, which is a tool to     %
  %  %   evaluate hemispherical images to estimate local sky view fraction  %
  %  %   and below-canopy shortwave radiation if the image is taken inside  %
  %  %   a forest. For more information see > help HPEval                   %
  %  %                                                                      %
  %  % USER SPECIFIC NOTES                                                  %
  %  %   This settings file is the default version for HPEval. Created by   %
  %  %   T.Jonas @ SLF/WSL, last modified on 2019/02/12                     %
  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% version of settings file
  % current version of settings file, any changes of the output structure 
  % array 'hpeset' (adding / removing a field, or changing the size / type
  % of an existing field) require creating a new version number so that 
  % necessary updates can be carried out upon loading hpeset within HPEval.
  % The versioning system of HPEval supports up to 3 levels, i.e. 1.4.2
  hpeset.version           = '1.2';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% access to input data
  % provide all paths/files in absolute terms 
  
  basefolder               = 'D:\Files\Matlab Software\ForestSnow\HPEval\';
  % use a basefolder to simplify the below path definitions
  
  hpeset.in.img            = fullfile(basefolder,'INPUT_IMAGES_BINARY\');
  % path to folder with hemispherical images to be analysed. Supported 
  % image formats include .jpg, .png, and .tif. For calculating below-
  % canopy shortwave radation, HPEval supports the use of a reference site 
  % that ideally represents above canopy conditions. The respective image 
  % file should be named reference*.*
  
  hpeset.in.inf            = fullfile(basefolder,'INPUT_IMAGES_BINARY\info.txt');
  % path to info file that contains specification of each image to be 
  % analysed including camera position, grayscale method, and evaluation 
  % period in case swr data is to be calculated. The format of the info 
  % file is available from the header of the example file. Note that time 
  % stamps must be given in the format dd.mm.yyyy HH:MM:SS, i.e. including
  % the second
  
  hpeset.in.dem            = fullfile(basefolder,'INPUT_DEM\dem_50m_crp.txt');
  % path to DEM data file. Leave empty to skip reading DEM, i.e.  
  % hpeset.in.dem = ''; Note however that a DEM is needed to a) mask 
  % surrounding terrain in the hemispherical image, and b) to evaluate 
  % local slope and aspect for the calculation of sub-canopy swr per 
  % inclined surface area
  
  hpeset.in.swr            = fullfile(basefolder,'INPUT_SWR_DATA');
  % path to folder with SWR data. Explicit files are to be listed in the 
  % info file. Note that SWR data from a reference site is needed to 
  % calculate true sub-canopy swr, rather than potential sub-canopy swr. In
  % case in.img / in.inf specify a refence site and associated swr data, 
  % the swr data is processed in two ways: a) any swr data corresponding to
  % a position of the sun behind canopy / terrain in the reference image is
  % removed, and b) the diffuse component of the SWR data is corrected 
  % should sky view fraction be < 1 at the reference site.
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% access to output data
  % provide all paths/files in absolute terms 
  
  hpeset.out.img           = fullfile(basefolder,'OUTPUT_IMAGES\');
  % path to folder for saving reclassified images. A respective folder is 
  % created should it not exist.
  
  hpeset.out.fig           = fullfile(basefolder,'OUTPUT_FIGURES\');
  % path to folder for saving resulting figures. A respective folder is 
  % created should it not exist.
  
  hpeset.out.dat           = fullfile(basefolder,'OUTPUT_DATA\');
  % path to folder for saving resulting swr data. A respective folder is 
  % created should it not exist.
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% define format of output data
  
  hpeset.for.img           = {'png'};                                       
  % specify output format of reclassified images, available options are:
  % 'png','tif','jpg'. Use {} to disable image output, multiple
  % entries are allowed. Put even a single entry into a cell array.
  
  hpeset.for.swr           = {'png'};
  % specify output format of swr figures, available options are: 
  % 'png','tif','jpg','fig'. Use {} to disable figure output, multiple
  % entries are allowed. Put even a single entry into a cell array.
  
  hpeset.for.dat           = {'mat'};
  % specify output format of swr data, available options are: 'mat','txt'.
  % Use {} to disable figure output, multiple entries are allowed. Put even
  % a single entry into a cell array.
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculation switches
  
  hpeset.sw.terrain        = 4;
  % this switch enables/disables masking of surrounding terrain. The 
  % following values are supported
  % 0 > will disabled masking, use this option in particular if no DEM is 
  %     available. This option may further speed up calculations. Note 
  %     however that enabling this option is helpful to avoid
  %     misclassification of snow on the ground as sky
  % 1 > to enable masking for every individual image
  % 2 > to enable masking for every individual image, but to speed up
  %     calculation by suppression the resampling of the horizon line to 
  %     an angular resolution of 1° (if calculated at coarser resolution).
  % 3 > to enable masking for the first image from which the resulting mask
  %     is inherited to any subsequent image. This option is recommended, 
  %     if the locations at which images were takes are, in relation to the 
  %     surrounding terrain, close to one another; and if images do not 
  %     feature considerably different rotational misalignments; and if the
  %     user has calculation time or CPU constraints. Note that this option
  %     entails the correction of rotational misalignments with regards to 
  %     the terrain mask also being inherited from the first image.
  % 4 > same as option 3, but plots terrain on top of canopy pixels, rather
  %     then below canopy pixels. This enables a much quicker transfer from
  %     the first to all subsequent images. While this option entaile an 
  %     unrealistic visual representation of terrain in the images, it has
  %     no implications for the resulting swr / svf values.
  % 5 > same as option 3, but uses terrain mask from first image excluding
  %     the reference image
  % 6 > same as option 4, but uses terrain mask from first image excluding
  %     the reference image
  
  hpeset.sw.svf            = 1;
  % this switch enables/disables calculation of sky view fraction. The 
  % following values are supported
  % 0 > will disable calculation of svf. Note that sky view fraction is 
  %     needed for the calculating of sub-canopy swr. If hpeset.sw.swr = 1
  %     svf will be calculated irrespective of this option, but in this 
  %     case no output file will be generated for svf data
  % 1 > will enable calculation of svf
  
  hpeset.sw.swr            = 1;
  % this switch enables/disables calculation of below-canopy shortwave
  % radation. The following values are supported
  % 0 > will disable calculation of swr
  % 1 > will enable calculation of swr, this is the recommended choice.
  % 2 > will enable calculation of swr, but will calculate the position of
  %     the sun only for the location and evaluation times of the first
  %     image (as specified in the info file) and inherit the sun track to 
  %     all subsequent images.
  
  hpeset.sw.detsun         = 0;
  % this switch enables/disables a detailed projection of the sun disk for 
  % the calculation of the transmissivity of the canopy for direct 
  % radiation. The following values are supported
  % 0 > will disable a detailed projection. In this case the projection of 
  %     the sun disk is approximated to take the form of a square. This
  %     will speed up the calculations massively with little impact on the
  %     resulting swr data. Therefor this option is highly recommended.
  % 1 > will enable a detailed projection. In this case the position, size
  %     and distortion of the projected sun disk is explcitly calculated
  %     for each image and time step. Unless maximum accuracy is needed,
  %     this option is not recommended.
  
  hpeset.sw.scatter        = 1;
  % this switch enables/disables simulating the effect of scattered direct 
  % swr due to haze, multiple reflections or refration within the canopy. 
  % The following values are supported
  % 0 > will disable simulating scattered direct swr   
  % 1 > will enable simulating scattered direct swr. Note that this option
  %     will entail increased calucation times
  
  hpeset.sw.suntrack       = 1;
  % this switch enables/disables inclusion of the track on sun disk in the 
  % output images. The following values are supported
  % 0 > will disable the inclusion of the sun track in the output images 
  % 1 > will enable the inclusion of the sun track in the output images 
  
  hpeset.sw.report         = 1;
  % this switch enables/disables the generation of command line messages 
  % after each processed image. The following values are supported
  % 0 > will disable the generation of command line messages 
  % 1 > will enable the generation of command line messages 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculation parameters

  hpeset.par.time_zone     = 1;
  % time zone used in data / info files, in offset hours relative to UTC, 
  % e.g. set to +1 if providing swr data from a data logger that uses time
  % stamps in CET.
  
  hpeset.par.coor_system   = 'CH1903';
  % coordinate system used in data / info files. While it is primarily 
  % important that all coordinates specified (DEM, camera position) feature
  % the same coordinate system, these need to be converted into lat/lon 
  % (WGS84) if calculation of swr is enabled. Only for this purpose the 
  % coordinate system needs to be specified. Supported option currently 
  % include
  % CH1903  > Swissgrid LV03
  % CH1903+ > Swissgrid LV95
  % UTM     > please note that in this case the utm zone must be specified
  
  hpeset.par.utm_zone      = '32N';
  % zone of utm coordinate system, using either the letter N or S as a 
  % specifier, whether or not the coordinates are given for the northern
  % or southern hemisfer, i.e. the letter shall not denote latitude bands
  % also sometimes used for specifying UTM coordinates
  
  hpeset.par.terrain_perim = 50000;
  % this parameter defines a horizontal perimeter around the camera
  % position within which the DEM is considered when calculating the
  % terrain mask. The perimeter is to be given in the same units as
  % the DEM (e.g meters in the case of UTM coordinates)
  
  hpeset.par.terrain_angs  = 90;
  % this parameter defines the number of azimuth angles for which the DEM
  % is evaluated to derive a horizon line for the terrain mask.
  
  hpeset.par.zenith_rings  = 9;
  % this parameter defines the number of zenith rings to be used when 
  % evaluating the hemispherical images for sky-view fraction.
  
  hpeset.par.lens_profile  = 1;
  % ID of a lens profile available in the database stored within HPEval, 
  % currently available profiles are: 
  % 0 > assuming a perfect equi-angular lens 
  % 1 > Sigma 4.5mm F2.8 EX DC HSM Circular Fisheye + Canon EOS 600D 
  % 2 > Yasuhara Madoka 180 7.3mm F4 Circular Fisheye + Sony Nex-6 (lens 1)
  % 3 > Yasuhara Madoka 180 7.3mm F4 Circular Fisheye + Sony Nex-6 (lens 2)
  % it is strongly recommended to create your own profile if working with
  % images from a different lens/camera combination
  
  hpeset.par.scatter_drad  = [0.533 1.066 2.132]./2;
  % this parameter defines the apparent size of the sun disk in degrees 
  % of the radius. While the true apparent diameter of the sun is about 
  % 0.533°, this number can be modified to conceptually account for 
  % scattered but yet directional swr (see hpeset.sw.scatter above). Up to 
  % three evaluation disks are supported by HPEval, so that the parameter 
  % can be a 1x1, 1x2, or 1x3 vector.
  
  hpeset.par.scatter_dwgh  = [50 30 20];
  % this parameter defines the relative weight each of the above sub disk
  % evaluation will be attributed to. The transmissivity of direct
  % radiation throught the canopy evaluated based on the first sun disk
  % scatter_drad(1) will be weighted with scatter_dwgh(1), and so on, to 
  % arrive at a weighted mean of all evaluation. Note that the the vector
  % scatter_dwgh must have at least the size of the vector scatter_drad.
  
  hpeset.par.swr_maxdir    = 3;
  % this parameter selects a particular way of constraining direct
  % radiation from the sun. A limiter is necessary to avoid unrealistic
  % radiation values in certain situations, such as calculating swr per
  % inclined surface area when the elevation angle of the sun is very low
  % and the terrain is facing towards the sun. Available choices are
  % 0 > using the solar constant
  % 1 > using measured swr transmissivity and modelled split into direct 
  %     and diffuse component
  % 2 > using de beer's law to acount for length of pathway in terms of 
  %     zenith angle. The extinction coefficient is based on the maximum 
  %     transmission of direct swr if sun is at zenith
  % 3 > applying constratins as in 1 and 2, which ever is stricter. This 
  %     is the default option

end