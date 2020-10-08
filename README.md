# strips_radiation_damage_2016  

This code makes the x-y distributions in the RP turning off U-V pairs. The code is based in the codes:  
https://github.com/cms-sw/cmssw/blob/master/RecoPPS/Local/plugins/TotemRPUVPatternFinder.cc  
https://github.com/cms-sw/cmssw/blob/master/RecoPPS/Local/plugins/TotemRPLocalTrackFitter.cc  

To run the code just folow the commands:

cmsrel CMSSW_10_6_12  
cd CMSSW_10_6_12/src/   
cmsenv  
git clone git@github.com:mrcthiel/strips_radiation_damage_2016.git  
scram b -j8  
cd strips_radiation_damage_2016/raddam_2016/  
cmsRun config_cfg.py name=fill5043  

This code is set to run over the skimmed data in:  
/eos/cms/store/group/phys_pps/reconstruction/2016/physics_runs/version-UL-2/   
the "name=fill5043" sets the file that will be used.  

In the output.root file there are the x-y distributions excluding a U-V pair. For example, the histogram x_y_0 is the x-y distribution excluding the planes 0 and 1, the histogram x_y_1 is the x-y distribution excluding the planes 2 and 3, the histogram x_y_2 is the x-y distribution excluding the planes 4 and 5... Histograms x_y_ > 4 are not used now.   

To count an event in the histogram x_y_0, the planes 0 and 1 (the ones that are turned off) need to have a consistent hit compared with the other planes. It is made in lines 613-766 in plugins/raddam_2016.cc. For example the hit in 0 needs at least match with one hit in 2 or 4 or 6 or 8, the U planes. To check the match, histograms doing all de differences (plane_0 - plane_2, plane_1 - plane_3, ...) are build previously (raddam_2016/parameters.root). The same is made for x_y_1, x_y_2, x_y_3 and x_y_4.  

# Running on condor  
An example to run this code using condor is in condor_sub/ directory. The to run over the 2016 DoubleEG C, the ListOfFiles_C.txt can be used. Before running condor, just compress the CMSSW and strips_radiation_damage_2016 directories in the files CMSSW.tar.gz and raddam_2016.tar.gz respectively.  


