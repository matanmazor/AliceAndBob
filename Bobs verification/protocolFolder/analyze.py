#analysis pipeline
#last updated on May 8th, 2017

from lib.preRNG import preRNG
import os
import random
import nibabel as nb


fsl_dir = "/usr/share/fsl"
template_dir = os.path.join(fsl_dir, "data", "standard")
atlas_dir = os.path.join(fsl_dir, "data", "atlases")
home_dir = os.getcwd()
stats_dir = os.path.join(home_dir, "analyzed", "all_runs.gfeat",
                             "cope1.feat", "stats")


def createEVs(protocol_folder):
    #set seed for pseudorandom number generator using the preRNG method
    preRNG(protocol_folder)
    os.mkdir(os.path.join(home_dir, "EVs")) #create EV directory
    for run_idx in [1,2,3,4]:        
        block_order = [0]*10+[1]*10
        random.shuffle(block_order)
        left_hand_EV = open(os.path.join(home_dir,"EVs", 
                                          "run"+str(run_idx)+"_LH"), "w")
        right_hand_EV = open(os.path.join(home_dir,"EVs", 
                                          "run"+str(run_idx)+"_RH"), "w")
        EVs = [left_hand_EV, right_hand_EV]
        time = 10 # run starts after ten seconds of rest
        for block in block_order:
            EVs[block].write("%d %d %d\n" % (time, 8, 1))
            time += 18 # each block+rest cycle is 18 seconds
        left_hand_EV.close()
        right_hand_EV.close()


def analyzeSingleRun(run_num):
    """
    distribute fsf files in functional dirs and run first level feat
    """
    infile = open(os.path.join(home_dir, "single_run.fsf"), "r")
    template = infile.read() 
    infile.close() 
    template = template.replace("run_num", str(run_num)) 
    template = template.replace("cur_dir", home_dir) 
    template = template.replace("templates_dir", template_dir)
    outfile = open(os.path.join(home_dir, "run"+str(run_num)+".fsf"),"w")
    outfile.write(template)
    outfile.close()
    
    os.system("fsl5.0-feat " + os.path.join(home_dir, 
                                "run" + str(run_num) + ".fsf"))
    

def analyzeAllRuns():
    """
    run feat on all four runs together (fixed effects analysis)
    """
    infile = open(os.path.join(home_dir, "all_runs_template.fsf"), "r") 
    template = infile.read() 
    infile.close() 
    template = template.replace("cur_dir", home_dir)
    template = template.replace("templates_dir", template_dir)
    outfile = open(os.path.join(home_dir, "all_runs.fsf") ,"w") 
    outfile.write(template)
    outfile.close()

    os.system("fsl5.0-feat " + os.path.join(home_dir, "all_runs.fsf"))



def FDR():
    """
    perform (two-sided) FDR correction for the cerebellum voxels
    following the steps described in fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDR
    """
    os.system(" ".join(["cp", "mask.nii.gz", stats_dir]))
    os.chdir(stats_dir)
    
    # create a p-value image

    dof = nb.load(os.path.join(stats_dir,"tdof_t1.nii.gz"))
    dof = dof.get_data().max()
    
    os.system(" ".join(["fsl5.0-ttologp -logpout logp1 varcope1 cope1",
                        str(dof)]))
    os.system("fsl5.0-fslmaths logp1 -exp p1")

    # perform fdr correction    
    pos_threshold = os.popen("fsl5.0-fdr -i p1 -m mask -q 0.025").read()
    pos_threshold = float(pos_threshold.split('\n')[1])
    os.system(" ".join(["fsl5.0-fslmaths p1 -mul -1 -add 1 -thr ",
                       str(1-pos_threshold),
                        " -mas mask thresh_1_minus_p1_pos"]))
    #apply mask to zstat1
    os.system("fsl5.0-fslmaths zstat1 -mas thresh_1_minus_p1_pos pos_fdr_zstat1")
    
    neg_threshold = os.popen("fsl5.0-fdr -i p1 -m mask -q 0.025 \
                                                    --oneminusp").read()
    neg_threshold = float(neg_threshold.split('\n')[1])
    os.system(" ".join(["fsl5.0-fslmaths p1 -thr ",
                       str(1-neg_threshold),
                        " -mas mask thresh_p1_neg"]))    
    #apply mask to zstat1
    os.system("fsl5.0-fslmaths zstat1 -mul -1 -mas thresh_p1_neg neg_fdr_zstat1")

    os.chdir(home_dir)


#create results directory
os.mkdir(os.path.join(home_dir, "analyzed"));

#create EV files
createEVs("path to protocol folder") #insert path to zipped protocol folder here
  
# SINGLE RUN
for run_num in [1,2,3,4]:
    analyzeSingleRun(run_num)

# ALL RUNS
analyzeAllRuns()

# SMALL_VOLUME FDR CORRECTION
FDR()

#visualize
os.system(" ".join(["fslview",
                    os.path.join(home_dir,"analyzed","run1.feat","reg",
                                     "highres2standard"),
                    os.path.join(stats_dir,"pos_fdr_zstat1"),
                    os.path.join(stats_dir,"neg_fdr_zstat1")]))

    
