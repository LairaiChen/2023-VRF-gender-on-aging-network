import os
import sys

#print("****************Entering New Job********************")
PATH_JOB = sys.argv[1] + "/"  # Gets working directory from command ars
os.chdir(PATH_JOB)
#print("Current Python Working Directory: " + PATH_JOB)

PATH_ATLAS_BN_246 = "atlases/BN_Atlas_246_1mm.nii"

PATH_ZIP_SOURCE = "biobank/images_all/"
PATH_RESULT = "coding/pipeline/camino/result/"

PATH_EID_TXT = PATH_JOB + "list.txt"


for eid in open(PATH_EID_TXT):
    os.chdir(PATH_JOB)
    eid = eid.strip()  # Gets rid of blank

    path_eid = PATH_JOB + eid + "/"

    path_image_t1_brain = path_eid + "t1/T1/T1_brain.nii.gz"
    path_image_nodif_brain = path_eid + "dti/dMRI/dMRI.bedpostX/nodif_brain.nii.gz"
    path_image_dti_to_t1 = path_eid  + "dti_to_t1.nii.gz"

    path_image_wmmask_t1 = path_eid + "t1/T1/T1_fast/T1_brain_pve_2.nii.gz"
    path_image_wmmask_dti = PATH_RESULT + "wm_mask_dti/" + eid + ".nii.gz"
    
    path_aal_to_t1 = path_eid + "aal_to_t1.nii.gz"
    path_aal_to_dti = PATH_RESULT + "atlas_aal_dti/" + eid + ".nii.gz"
    path_bn_to_t1 = path_eid + "bn_to_t1.nii.gz"
    path_bn_to_dti = PATH_RESULT + "atlas_bn_dti/" + eid + ".nii.gz"

    path_mat_dti_to_t1 = path_eid + "dti_to_t1.mat"
    path_mat_t1_to_dti = path_eid + "t1_to_dti.mat"

    path_warp_t1_to_mni = path_eid + "t1/T1/transforms/T1_to_MNI_warp_coef.nii.gz"
    path_warp_mni_to_t1 = path_eid + "mni_to_t1_warp.nii.gz"

    path_camino_track = path_eid + eid + "_camino_bedpost_track"
    path_camino_track_post = PATH_RESULT + "tracks_post/" + eid

    '''1. Makes working dir and unzips'''
    #start = time.time()
    command_mkdir = "mkdir" + " " + eid
    os.system(command_mkdir)

    #print("Unzipping dti for " + eid)
    # Unzip dti
    command_unzip_dti = "unzip -q " + PATH_ZIP_SOURCE + eid + "_20250_2_0.zip"  + " -d "  +  path_eid + "dti/"
    
    # Unzip T1WI
    #print("Unzipping t1 for " + eid)
    command_unzip_t1 = "unzip -q " + PATH_ZIP_SOURCE + eid + "_20252_2_0.zip" + " -d "  +  path_eid+ "t1/"
    # print(command_unzip_dti)
    os.system(command_unzip_dti)

    # print(command_unzip_t1)
    os.system(command_unzip_t1)

    #end = time.time()
    #print("Unzipping takes: " + str(end-start) + 's')

    '''2. Registers atlas and WM mask'''
    #start = time.time()

    '''2.1 Creates transformations'''
    #print("Creating transformations for " + eid)
    # Creates warp fields from mni to t1
    command_invwap_mni_to_t1 = "invwarp --ref=" + path_image_t1_brain + " --warp=" + path_warp_t1_to_mni + " --out=" + path_warp_mni_to_t1
    #print(command_invwap_mni_to_t1)
    os.system(command_invwap_mni_to_t1)

    # Creates transformation mat from t1 to dti

    # First register dti to t1
    command_flirt_mat_dti_to_t1 = "flirt" + " -in " + path_image_nodif_brain + " -ref " + path_image_t1_brain + " -out " + path_image_dti_to_t1 + " -omat " + path_mat_dti_to_t1 + " -interp nearestneighbour"
     # Then inverse the mat
    command_inverse_mat_dti_to_t1 = "convert_xfm" + " -omat " + path_mat_t1_to_dti + " -inverse " + path_mat_dti_to_t1
    #print(command_flirt_mat_dti_to_t1)
    #print(command_inverse_mat_dti_to_t1)
    os.system(command_flirt_mat_dti_to_t1)
    os.system(command_inverse_mat_dti_to_t1) 

    '''2.2 applies transformations'''
    # Register bn246 (mni -> t1 -> dti)
    #print("Regestering BN246 for " + eid)
    path_atlas = PATH_ATLAS_BN_246
    path_atlas_to_t1 = path_bn_to_t1
    path_atlas_to_dti = path_bn_to_dti
    command_atlas_mni_to_ti = "applywarp" + " -i " + path_atlas + " -r " +  path_image_t1_brain + " -w " + path_warp_mni_to_t1 + " -o " + path_atlas_to_t1 + " --interp=nn"
    command_atlas_t1_to_dti = "flirt" + " -in " + path_atlas_to_t1 + " -init " + path_mat_t1_to_dti + " -ref " + path_image_nodif_brain+ " -out " + path_atlas_to_dti + " -applyxfm" + " -interp nearestneighbour"
    
    #print(command_atlas_mni_to_ti)
    os.system(command_atlas_mni_to_ti)
    #print(command_atlas_t1_to_dti)
    os.system(command_atlas_t1_to_dti)

    #print(command_atlas_mni_to_ti)
    os.system(command_atlas_mni_to_ti)
    #print(command_atlas_t1_to_dti)
    os.system(command_atlas_t1_to_dti)
    
    # Register WM mask
    command_wmmask_from_t1_to_dti = "flirt" + " -in " + path_image_wmmask_t1 + " -init " + path_mat_t1_to_dti + " -ref " + path_image_nodif_brain+ " -out " + path_image_wmmask_dti + " -applyxfm" + " -interp nearestneighbour"
    #print(command_wmmask_from_t1_to_dti)
    os.system(command_wmmask_from_t1_to_dti)
    #end = time.time()
    #print("Registration takes: " + str(end-start) + 's')


    '''3. Run CAMINO'''
    #start = time.time()
 
    # Track with bedpost and rk4
    command_camino_track = "track -inputmodel bedpostx_dyad -curvethresh 45 -curveinterval 5 -bedpostxminf 0.1  -header " + path_image_nodif_brain + " -seedfile " + path_image_wmmask_dti + " -bedpostxdir " + path_eid + "dti/dMRI/dMRI.bedpostX -tracker rk4 -interpolator nn -stepsize 2 -outputfile " + path_camino_track
    #print(command_camino_track)
    os.system(command_camino_track)

    # Post process streamlines
    command_camino_track_post = "procstreamlines -inputfile " + path_camino_track + " -mintractlength 20 -maxtractlength 250 -header " + path_image_nodif_brain + "  > " + path_camino_track_post
    #print(command_camino_track_post)
    os.system(command_camino_track_post)

    # Construct SC-network based on BN246
    path_atlas_dti = path_bn_to_dti
    path_output_csv = path_eid + "csv_raw_bn_"
    command_sc_atlas = " conmat -inputfile " + path_camino_track_post + " -targetfile " + path_atlas_dti + " -outputroot " + path_output_csv
    os.system(command_sc_atlas)

    #end = time.time()
    #print("CAMINO takes: " + str(end-start) + 's')

    '''Generates csv and edge files'''
    # Remove header and set diag to zero

    '''BN'''      
    path_mat_raw = path_eid + "csv_raw_bn_" + 'sc.csv'
    path_mat_header_rm = PATH_RESULT + "sc_bn_header/" + eid + '.csv'
    path_mat_diag =  PATH_RESULT + "sc_bn_diag/" + eid + '.csv'
    command_rm_header = 'Rscript coding/transforms/csv/csv_remove_header.R ' + path_mat_raw + " " + path_mat_header_rm
    os.system(command_rm_header)
    command_diag_mat = 'Rscript coding/transforms/csv/csv_diag_zero.R ' + path_mat_header_rm + " " + path_mat_diag
    os.system(command_diag_mat)

    '''Saves the results'''
    # Masks, atalases, are saved earlier in this code
    os.system("cp " + path_image_t1_brain + " " + PATH_RESULT + "t1wi/" + eid + ".nii.gz") # T1WI
    os.system("cp " + path_eid + "dti/dMRI/dMRI/dti_FA.nii.gz" + " " + PATH_RESULT + "dti_fa/" + eid + ".nii.gz") # fa
    os.system("cp " + path_eid + "dti/dMRI/dMRI/dti_MD.nii.gz" + " " + PATH_RESULT + "dti_md/" + eid + ".nii.gz") # md
    os.system("cp " + path_eid + "dti/dMRI/dMRI/NODDI_ICVF.nii.gz" + " " + PATH_RESULT + "noddi_icvf/" + eid + ".nii.gz") # icvf
    os.system("cp " + path_eid + "dti/dMRI/dMRI/NODDI_ISOVF.nii.gz" + " " + PATH_RESULT + "noddi_isovf/" + eid + ".nii.gz") # isovf
    os.system("cp " + path_eid + "dti/dMRI/dMRI/NODDI_OD.nii.gz" + " " + PATH_RESULT + "noddi_od/" + eid + ".nii.gz") # od
    os.system("cp " + path_image_nodif_brain + " " + PATH_RESULT + "dti_b0/" + eid + ".nii.gz") # b0
    os.system("cp " + path_eid + "dti/dMRI/dMRI/bvecs" + " " + PATH_RESULT + "dti_bvecs/" + eid) # bvecs
    os.system("cp " + path_eid + "dti/dMRI/dMRI/bvals" + " " + PATH_RESULT + "dti_bvals/" + eid) # bvals
    os.system("cp " + path_warp_t1_to_mni + " " + PATH_RESULT + "warp_t1_to_mni/" + eid + ".nii.gz") # Warp field from t1 to mni space
    os.system("cp " + path_warp_mni_to_t1 + " " + PATH_RESULT + "warp_mni_to_t1/" + eid + ".nii.gz") # Warp field from mni to t1 space
    os.system("cp " + path_mat_dti_to_t1 + " " + PATH_RESULT + "mat_dti_to_t1/" + eid + ".mat") # Linear transformation mat from dti to t1
    os.system("cp " + path_mat_t1_to_dti + " " + PATH_RESULT + "mat_t1_to_dti/" + eid + ".mat") # Linear transformation mat from t1 to dti
    os.system("cp -r " + path_eid + "t1/T1/T1_fast" + " " + PATH_RESULT + "fast/" + eid) # fast
    os.system("cp -r " + path_eid + "t1/T1/T1_first" + " " + PATH_RESULT + "first/" + eid) # fast


    '''Free spaces'''
    os.chdir(path_eid)
    os.system("rm -r -f *")
    #print("************************Done for " + eid + "!********************************")
