#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 06:18:53 2017

@author: Bob

preRNG verification for the paper "Cerbellum Involvement in Hand Movements: 
    a functional Magnetic Resonance Imaging Experiment" by Alice, 2017. 
    Full documentation is available in the accompanying PDF file.
    
"""
from lib.preRNG import preRNG
import os
import random
import nibabel as nib
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import stats

fsl_dir = "/usr/share/fsl"
template_dir = os.path.join(fsl_dir, "data", "standard")
atlas_dir = os.path.join(fsl_dir, "data", "atlases")
home_dir = os.getcwd()
stats_dir = os.path.join(home_dir, "analyzed", "all_runs.gfeat",
                             "cope1.feat", "stats")



def changeProtocol():
    with open(os.path.join(home_dir,"protocolFolder","analyze.py"), "a") as f:
        f.write("#")
        f.close()
    os.remove(os.path.join(home_dir,"protocolFolder.zip"))
    shutil.make_archive(os.path.join(home_dir,"protocolFolder"),
                        'zip', os.path.join(home_dir,"protocolFolder"))

def createEVs(protocol_folder):
    shutil.rmtree(os.path.join(home_dir, "EVs"))
    #set seed for pseudorandom number generator using the preRNG method
    protocolSum = preRNG(protocol_folder)
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
    return protocolSum

def ppSingleRun(run_num):
    """
    distribute fsf files in functional dirs and run first level feat
    """
    infile = open(os.path.join(home_dir, "pp_single_run.fsf"), "r")
    template = infile.read() 
    infile.close() 
    template = template.replace("run_num", str(run_num)) 
    template = template.replace("cur_dir", home_dir) 
    template = template.replace("template_dir", template_dir)
    outfile = open(os.path.join(home_dir, "run"+str(run_num)+".fsf"),"w")
    outfile.write(template)
    outfile.close()
    
    os.system("fsl5.0-feat " + os.path.join(home_dir, 
                                "run" + str(run_num) + ".fsf"))

def analyzeSingleRun(run_num, iteration):
    """
    distribute fsf files in functional dirs and run first level feat
    """
    shutil.copytree(os.path.join(home_dir,'analyzed','run'+str(run_num)+'.feat'),
                    os.path.join(home_dir,'analyzedWrongPF',str(iteration),
                                 'run'+str(run_num)+'.feat'))
    infile = open(os.path.join(home_dir, "single_run.fsf"), "r")
    template = infile.read() 
    infile.close() 
    template = template.replace("run_num", str(run_num)) 
    template = template.replace("iteration", str(iteration)) 
    template = template.replace("cur_dir", home_dir) 
    template = template.replace("template_dir", template_dir)
    outfile = open(os.path.join(home_dir, "run"+str(run_num)+".fsf"),"w")
    outfile.write(template)
    outfile.close()
    
    os.system("fsl5.0-feat " + os.path.join(home_dir, 
                                "run" + str(run_num) + ".fsf"))


def analyzeAllRuns(iteration):
    """
    run feat on all four runs together (fixed effects analysis)
    """
    infile = open(os.path.join(home_dir, "all_runs_template.fsf"), "r") 
    template = infile.read() 
    infile.close() 
    template = template.replace("cur_dir", home_dir)
    template = template.replace("iteration", str(iteration)) 
    template = template.replace("template_dir", template_dir)
    outfile = open(os.path.join(home_dir, "all_runs.fsf") ,"w") 
    outfile.write(template)
    outfile.close()

    os.system("fsl5.0-feat " + os.path.join(home_dir, "all_runs.fsf"))

def takeSnapshots(iteration):        
    os.system(" ".join(["fsl5.0-slicer", 
        os.path.join(home_dir,"analyzedWrongPF",str(iteration),"all_runs.gfeat",
                         "cope1.feat", "rendered_thresh_zstat1.nii.gz"),
        "-z 0.7",
        os.path.join(home_dir,"analyzedWrongPF",
                        str(iteration), "axial"),
        "-x 0.7",
        os.path.join(home_dir,"analyzedWrongPF",
                        str(iteration), "saggital"),
        "-y 0.45",
        os.path.join(home_dir,"analyzedWrongPF",
                        str(iteration), "coronal")]))
        
def computeLikelihood(iteration,volumes_per_run):
    #compute pooled variance
    variance = 0
    sum_squared = 0
    for run in [1,2,3,4]:
        variance+=nib.load(os.path.join(home_dir,"analyzedWrongPF",str(iteration),
                                    "run"+str(run)+".feat", "stats",
                                    "sigmasquareds.nii.gz")).get_data()
        sum_squared+= np.power(nib.load(os.path.join
                                        (home_dir,"analyzedWrongPF",str(iteration),
                                         "run"+str(run)+".feat", "stats",
                                         "res4d.nii.gz")).get_data(),
                                2).sum(3)
    mask = variance>0
    variance/=4
    loglikelihood = np.zeros(variance.shape)
    loglikelihood[mask] = -(volumes_per_run*4)/2*np.log(variance[mask])-sum_squared[mask]/(2*variance[mask])

def createFunctionalMasksFromFunc():
    
    example_func = nib.load(os.path.join("analyzed","run1.feat","reg","example_func.nii.gz"))
    
    
    left_S1_centroid = [68,37,54] #in functional space, voxel coordinates
    leftS1 = np.zeros(example_func.shape, 'int32')
    leftS1[left_S1_centroid[0]-1:left_S1_centroid[0]+2,
           left_S1_centroid[1]-1:left_S1_centroid[1]+2,
           left_S1_centroid[2]-1:left_S1_centroid[2]+2]=1
    img = nib.Nifti1Image(leftS1, example_func.affine)
    img.to_filename('leftS1.nii.gz')
    
    right_S1_centroid = [30, 37, 54]
    rightS1 = np.zeros(example_func.shape, 'int32')
    rightS1[right_S1_centroid[0]-1:right_S1_centroid[0]+2,
            right_S1_centroid[1]-1:right_S1_centroid[1]+2,
            right_S1_centroid[2]-1:right_S1_centroid[2]+2]=1
    img = nib.Nifti1Image(rightS1, example_func.affine)
    img.to_filename('rightS1.nii.gz')
    
    ## transform to MNI space
    os.system(" ".join(["fsl5.0-applywarp",
                        "--ref="+os.path.join("analyzed","run1.feat","reg","standard"),
                        "--in="+"leftS1",
                        "--warp="+os.path.join("analyzed","run1.feat","reg",
                                               "highres2standard_warp"),
                        "--premat="+os.path.join("analyzed","run1.feat","reg",
                                                  "example_func2highres.mat"),
                        "--out="+"leftS1standard.nii"]))
        
    os.system(" ".join(["fsl5.0-applywarp",
                        "--ref="+os.path.join("analyzed","run1.feat","reg","standard"),
                        "--in="+"rightS1",
                        "--warp="+os.path.join("analyzed","run1.feat","reg",
                                               "highres2standard_warp"),
                        "--premat="+os.path.join("analyzed","run1.feat","reg",
                                                  "example_func2highres.mat"),
                        "--out="+"rightS1standard.nii"]))
    
def computeLogLikelihood(right_vec, left_vec, right_timecourse, left_timecourse):
    #create true regressor files
    (right_regressor, left_regressor) = (np.zeros((185,)), np.zeros((185,)))
    for block in left_vec:
         left_regressor[block/2+1:block/2+5]=1
    for block in right_vec:
        right_regressor[block/2+1:block/2+5]=1    
    right_hand_model = np.polyfit(right_regressor,left_timecourse,1,full=True)
    left_hand_model = np.polyfit(left_regressor,right_timecourse,1,full=True)
    right_sigma_squared = right_hand_model[1]/185
    left_sigma_squared = left_hand_model[1]/185
    true_log_likelihood = -(185/2)*(np.log(4*np.pi)
                            +np.log(left_sigma_squared)
                            +np.log(right_sigma_squared)
                            +2)   
    return true_log_likelihood

def swapBlocks(right_vec, left_vec, i, j):
    (swapped_right_vec, swapped_left_vec) = (right_vec[:],left_vec[:])
    swapped_right_vec.append(swapped_left_vec.pop(i))
    swapped_left_vec.append(swapped_right_vec.pop(j))
    swapped_right_vec.sort()
    swapped_left_vec.sort()
    return (swapped_right_vec, swapped_left_vec)

def gradientDescent(right_timecourse, left_timecourse):
    #initialize right and left vecs
    order_vec = range(10,370,18)
    random.shuffle(order_vec)
    (right_vec, left_vec) = (order_vec[:10], order_vec[10:])
    right_vec_mat = []
    myiter = 0
    max_iter = 1000
    fig = plt.figure(figsize = (25,5))
    ax1 = plt.subplot2grid((1, 6), (0, 0), colspan=5)
    ax1.set_ylim([-2.5, 5.5])
    ax1.plot(left_timecourse, linewidth = 3, color = "#0178B3") #blue
    ax1.plot(right_timecourse, linewidth = 3, color = "#D13F35") #red
    ax1.set_yticks([])
    model_log_likelihood = computeLogLikelihood(right_vec, left_vec, 
                                                right_timecourse, left_timecourse)
    found_local_maximum = False
    while not found_local_maximum and myiter<max_iter:
        BF_mat = obtainBFmat(right_vec, left_vec, right_timecourse, left_timecourse)
        model_log_likelihood = computeLogLikelihood(right_vec, left_vec, 
                                                right_timecourse, left_timecourse)
        ax2 = plt.subplot2grid((1, 6), (0, 5))
        plotBFmat(BF_mat, ax2)
        ax1.text(-3, 5-(myiter*0.7),"%.0f" 
                 % model_log_likelihood[0], fontsize =17)
        for block in left_vec:
            ax1.add_patch(patches.Rectangle(
                    (block/2, 5-(myiter*0.7)),   # (x,y)
                    4,          # width
                    0.3,          # height
                    color = "#D13F35"))
        for block in right_vec:
            ax1.add_patch(patches.Rectangle(
                    (block/2, 5-(myiter*0.7)),   # (x,y)
                    4,          # width
                    0.3,          # height
                    color = "#0178B3"))
        if myiter%100 == 0:
            print("for iteration number %d the highest logBF was %f" %(myiter, np.max(BF_mat)))
        found_local_maximum = np.max(BF_mat)<=0                            
        best_proposal = np.where(BF_mat == np.max(BF_mat))
        if not found_local_maximum:
            #draw parallelograms
            print(left_vec)
            
            left_swap = left_vec.pop(best_proposal[0])
            right_swap = right_vec.pop(best_proposal[1])
            print best_proposal[0], best_proposal[1]            
            ax1.add_patch(patches.Rectangle(
                    (left_swap/2, 5-(myiter*0.7)),   # (x,y)
                    4,          # width
                    0.3,          # height
                    fill = None, alpha = 1, linewidth = 3))
            ax1.add_patch(patches.Rectangle(
                    (right_swap/2, 5-(myiter*0.7)),   # (x,y)
                    4,          # width
                    0.3,          # height
                   fill = None, alpha = 1, linewidth = 3))
            ax2.scatter(best_proposal[1], best_proposal[0],
                        marker='x',color='white')
            right_vec.append(left_swap)
            left_vec.append(right_swap)
            right_vec_mat.append(right_vec[:])
        myiter+=1
        fig.savefig("_".join(["GD",str(run_num),str(myiter)])+".png", dpi=300)
        fig.delaxes(ax2)
    return (right_vec, right_vec_mat)

def obtainBFmat(right_vec, left_vec, right_timecourse, left_timecourse):
    right_vec.sort()
    left_vec.sort()
    BF_mat = np.zeros((10,10))
    true_log_likelihood = computeLogLikelihood(right_vec,left_vec, 
                                               right_timecourse, left_timecourse) 
    for i in range(10):
        for j in range(10):
            (swpd_right_vec, swpd_left_vec) = swapBlocks(right_vec,left_vec,i,j)
            log_likelihood = computeLogLikelihood(swpd_right_vec, swpd_left_vec,
                                                  right_timecourse, left_timecourse)
            BF_mat[i,j] = log_likelihood-true_log_likelihood  
    return BF_mat

def plotBFmat(BF_mat, ax):
    max_abs_BF = np.max(np.abs(BF_mat))
    img = ax.imshow(BF_mat, clim = (-max_abs_BF,max_abs_BF), cmap="PiYG")
    plt.colorbar(img, ax = ax) 
    ax.set_xticks(np.arange(0, 10, 1));
    ax.set_yticks(np.arange(0, 10, 1));
    ax.set_xticklabels(np.arange(1, 11, 1));
    ax.set_yticklabels(np.arange(1, 11, 1));
    # Minor ticks
    ax.set_xticks(np.arange(-.5, 10, 1), minor=True);
    ax.set_yticks(np.arange(-.5, 10, 1), minor=True);

    #Gridlines based on minor ticks
    ax.grid(which='minor', color='w', linestyle='-', linewidth=2)

    ax.set_xlabel('right hand')
    ax.set_ylabel('left hand')
    plt.savefig("run"+str(run_num)+"BF.png", dpi=300)    






        
#### Steps 1-3 as described in the PDF file
for run_num in [1,2,3,4]:
    ppSingleRun(run_num)

os.mkdir(os.path.join(home_dir, "analyzedWrongPF")) 
os.mkdir(os.path.join(home_dir, "EVs"))
#make a backup copy of analyze.py before making changes
shutil.copy2(os.path.join(home_dir, "protocolFolder", "analyze.py"),
             os.path.join(home_dir, "original_analyze.py"))

for iteration in range(100):  
    print("iter %d" % iteration)
    os.mkdir(os.path.join(home_dir, "analyzedWrongPF", str(iteration))) 
    shutil.copy(os.path.join(home_dir,"protocolFolder","analyze.py"), 
                os.path.join(home_dir, "analyzedWrongPF", str(iteration),"analyze.py")) 
    protocolSum = createEVs(os.path.join(home_dir,"protocolFolder.zip"))
    f = open(os.path.join(home_dir, "analyzedWrongPF", str(iteration),"sum.txt"),"w")  
    f.write(protocolSum)
    f.close() 
    print("created directory and preRNGed")
    for run_num in [1,2,3,4]:
        analyzeSingleRun(run_num, iteration)  
        print("analyzed run %d" % run_num)
    analyzeAllRuns(iteration)
    print("analyzing all runs")
    takeSnapshots(iteration)
    print("taking snapshots")
    for run_num in [1,2,3,4]:
        shutil.rmtree(os.path.join(home_dir,"analyzedWrongPF",str(iteration), 
        "run"+str(run_num)+".feat"))
    print("deleted all feat dirs")
    changeProtocol()
    print("changed protocol")

# Steps 4-6 as described in the PDF file
#set the protocol folder back to its original state and create new EVs
shutil.copy2(os.path.join(home_dir, "original_analyze.py"),
             os.path.join(home_dir, "protocolFolder", "analyze.py"))
os.remove(os.path.join(home_dir,"protocolFolder.zip"))
shutil.make_archive(os.path.join(home_dir,"protocolFolder"),
                        'zip', os.path.join(home_dir,"protocolFolder"))
createEVs(os.path.join(home_dir,"protocolFolder.zip"))

#extract the timecourse of the right and left primary motor areas
createFunctionalMasksFromFunc()
os.system(" ".join(["fsl5.0-invwarp",
                        "--ref="+os.path.join("analyzed","run1.feat","reg","highres"),
                        "--warp="+os.path.join("analyzed","run1.feat","reg","highres2standard_warp"),
                        "--out="+"highres2standard_warp_inv"]))

for run_num in [1,2,3,4]:
    
    # transform masks to native functional space, based on the steps found in:
    # https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FNIRT/UserGuide#Transforming_an_image_.28.27.27e.g..27.27_a_mask.29_in_standard_space_into_functional_space
    os.system(" ".join(["fsl5.0-convert_xfm",
                        "-omat", "highres2example_func.mat",
                        "-inverse", os.path.join("analyzed","run"+str(run_num)+".feat","reg","example_func2highres.mat")]))
        
    os.system(" ".join(["fsl5.0-applywarp",
                            "--ref="+os.path.join("analyzed","run"+str(run_num)+".feat","reg","example_func"),
                            "--in="+"leftS1standard",
                            "--warp="+"highres2standard_warp_inv",
                            "--postmat=highres2example_func.mat",
                            "--out="+"leftS1func"]))
        
    os.system(" ".join(["fsl5.0-fslmaths",
                            "leftS1func",
                            "-thr", "0.5",
                            "-bin", "leftS1func"]))
        
    os.system("fsl5.0-fslchfiletype NIFTI leftS1func")
        
    os.system(" ".join(["fsl5.0-applywarp",
                            "--ref="+os.path.join("analyzed","run"+str(run_num)+".feat","reg","example_func"),
                            "--in="+"rightS1standard",
                            "--warp="+"highres2standard_warp_inv",
                            "--postmat=highres2example_func.mat",
                            "--out="+"rightS1func"]))
        
    os.system(" ".join(["fsl5.0-fslmaths",
                            "rightS1func",
                            "-thr", "0.5",
                            "-bin", "rightS1func"]))
        
    os.system("fsl5.0-fslchfiletype NIFTI rightS1func")
    os.system("fsl5.0-fslchfiletype NIFTI "+
              os.path.join("analyzed","run"+str(run_num)+".feat","filtered_func_data"))
    img = nib.load(os.path.join("analyzed","run"+str(run_num)+".feat","filtered_func_data.nii"))
    img = img.get_data(caching='unchanged')
    right_S1 = nib.load("rightS1func.nii")
    right_S1 = right_S1.get_data()
    right_S1_voxels = np.ma.array(img, mask=img*(1-right_S1)[:,:,:, np.newaxis])
    right_timecourse = stats.zscore(right_S1_voxels.mean(axis=0).mean(axis=0).mean(axis=0).data)
    left_S1 = nib.load("leftS1func.nii")
    left_S1 = left_S1.get_data()
    left_S1_voxels = np.ma.array(img, mask=img*(1-left_S1)[:,:,:, np.newaxis])
    left_timecourse = stats.zscore(left_S1_voxels.mean(axis=0).mean(axis=0).mean(axis=0).data)
    
    #obtain and plot EVs:
    fname = "EVs/run"+str(run_num)+"_LH"
    with open(fname,"r") as f:
        left_vec = [int(r.split()[0]) for r in f]
    fname = "EVs/run"+str(run_num)+"_RH"
    with open(fname,"r") as f:
        right_vec = [int(r.split()[0]) for r in f]
    
    fig = plt.figure(figsize = (40,10))
    ax1 = fig.add_subplot(111)
    ax1.set_yticks([])
    for block in left_vec:
        ax1.add_patch(patches.Rectangle(
        (block/2, 3),   # (x,y)
        4,          # width
        0.3,          # height
        color = "#D13F35"))
    for block in right_vec:
        ax1.add_patch(patches.Rectangle(
        (block/2, 3),   # (x,y)
        4,          # width
        0.3,          # height
        color = "#0178B3"))
    
    ax1.plot(left_timecourse, linewidth = 6, color = "#0178B3") #blue
    ax1.plot(right_timecourse, linewidth = 6, color = "#D13F35") #red
    plt.savefig("run"+str(run_num)+".png")

    #plot BF matrix
    fig = plt.figure(figsize = (8,8))
    ax1 = fig.add_subplot(111)
    BF_mat = obtainBFmat(right_vec, left_vec, right_timecourse, left_timecourse)
    plotBFmat(BF_mat, ax1)
    plt.savefig("BFmat_run"+str(run_num)+".png")

    gradientDescent(right_timecourse, left_timecourse)    
    
