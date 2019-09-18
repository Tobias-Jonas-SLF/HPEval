 ### GENERAL DESCRIPTION                                                  
   This scripts evaluates real hemispherical images to estimate local sky view fraction and below-canopy shortwave radiation if the image is taken inside a forest.                                    
                                                                      
 ### VERSION                                                              
   This is a first complete version of HPEval to be references as vs 1.2 | created 2018/08 | last modified on 2019/09/14             
                                                                      
 ### CONTRIBUTING AUTHORS                                                 
   T.Jonas<sup>1)</sup>, C.Webster<sup>1,2)</sup>, G.Mazzotti<sup>1)</sup>, J.Malle<sup>1,3)</sup>            
   <sup>1)</sup> WSL Institute for Snow and Avalanche Research SLF, Davos, CH  
   <sup>2)</sup> University of Edinburgh, School of GeoSciences, Edinburgh, UK   
   <sup>3)</sup> Northumbria University, Newcastle upon Tyne, UK               
                                                                      
 ### CITATION (to be updated)                                             
   T.Jonas, C.Webster, G.Mazzotti, and J.Malle; 2019; HPEval: A canopy shortwave transmission model using high resolution hemispherical images; Agricultural and Forest Meteorology, in review                                                          
                                                                      
 ### CONTRIBUTING MATERIAL                                                
   - solar position calculations from NOAA, found online at www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html (accessed 2018/06)                                               
   - function to split solar radiation into direct/diffuse component as in the factorial snow model FSM 1.0 by R. Essery available at www.geosci-model-dev.net/8/3867/2015 (accessed 2018/06)                                               
   - function to convert UTM into WGS84 coordinates, found online at www.mathworks.com/matlabcentral/fileexchange/44242-utm2lonlat (accessed 2018/08)                                               
                                                                      
 ### TEST DATASET                                                         
   - hemispherical images taken around Davos Switzerland and corresponding meta data; by SLF                                  
   - shortwave radiation data measured at locations were above images were taken test images; by SLF                                     
   - DEM at 100 m spatial resolution, derived from TopoSwiss data available for research; by TopoSwiss                             
                                                                      
 ### OUTPUT                                                               
   The script has a number of optional outputs, which include         
   - reclassification of image (canopy, terrain, sky, solar track)    
   - shortwave radiation (per horizontal, inclined, sun-facing area)  
   - sky view fraction (for flat and/or hemispherical sensor)         
                                                                      
 ### METHODS                                                              
   The scripts includes the following steps (some of them being optional):                                                         
   - grayscale image (according to external input on threshold and thresholding method)                                             
   - correct image distortion (using a lens-specific profile)         
     note: this is an essential feature and should not be skipped     
   - mask surrounding terrain and add to image (based on DEM)         
     note: helps to avoid misclassification of snow on the ground     
   - calculate sky view fraction                                      
   - calculate position of sun (according to formulas by NOAA)        
   - plot sun disk on image behind canopy / terrain elements to evaluate the transmissivity of direct radiation through the canopy                                                           
   - split swr data measured at open site into diffuse/direct radiation (according to formulas in FSM 1.0 by R. Essery)        
     note: potential swr is used if no measurements are available     
   - calculate sub-canopy swr based on diffuse radiation times sky view fraction plus direct radiation times transmissivity         
                                                                      
 ### NOTATION                                                             
   local spherical coordinates are given as [rad,tht,phi], where     
   - rad is radial distance in [m]                                    
   - tht is zenith angle (theta) in degree, i.e.                      
       sun @ zenith is tht = 0;                                       
       sun @ horizon is tht = 90                                      
   - phi is azimuthal angle (phi) in degree, where                    
       0° is North                                                    
       90° is East                                                    
       180° is South                                                  
       270° is West                                                   
       note that hemispherical images are taken upwards so that azimuthal angles increase in counterclockwise direction, i.e. South is up if West is left                               
                                                                      
 ### SETTINGS                                                             
   All settings such as paths to I/O data and model parameters are handled through an external settings file the path/filename of which is the only input argument to this script                    
                                                                      
 ### USAGE                                                                
   run from the Matlab command line: HPEval('HPEval_Settings.m')                                      
