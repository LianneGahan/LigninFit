import numpy as np
import sys

def findmin(x,y):
    if x <= y:
        return x
    else:
        return y


vardata = np.loadtxt("vars.txt")
test_shape = np.full((1,2),0.)

#print("VARDATA SHAPE: " + str(vardata.shape))
noise_range = 0.;
if vardata.size == 2:
    mini_posi = np.full(1,1)
    np.savetxt("best_candidate.txt", mini_posi, delimiter = "\t", fmt = "%1i")    
else:
    mini_nr = vardata[0,0]
    mini = vardata[0,1]
    oldmini = mini

    # Find the smallest variance in this generation
    for i in range(1,vardata[:,0].size):
        oldmini = mini
        mini = findmin(mini, vardata[i,1])
        if mini < oldmini:
            mini_nr = vardata[i,0]
    print("MINI NR IS: ", mini_nr)   
    print("MINI IS: ", mini)     

    previous_mini = np.loadtxt("smallest_var.txt")#This is the minimum var from the generation before
    count_unsuccessfull = np.loadtxt("no_reduction_count.txt")
    count_unsuccessfull_out = np.full(1,count_unsuccessfull)
    #print("Old min variance: " + str(previous_mini))

    if mini < previous_mini:
#        print ("Old mini = " + str(previous_mini) + "; new mini = " + str(mini))
        # Read in best kin specs, which now becomes the second best or "old best"
        old_kin_data = np.loadtxt("best_kin_specs.txt");
        with open("old_best_kin_specs.txt", "w") as f:
            for item in old_kin_data:
                f.write("%1.8f\t" % item)
        # Also write the second best variance 
        with open("old_smallest_var.txt", "w") as f:
            f.write("%1.8f\t" % previous_mini)


        new_kin_data = np.loadtxt("family_" + str(sys.argv[1]) + "/Generation_" + str(sys.argv[2]) + "/Run_" + str(int(mini_nr)) + "/Params/kinetic_parameters.txt")
        new_kin_data_out = np.full((1,new_kin_data[:].size),0.)
        new_kin_data_out[0,:] = new_kin_data[:]
        np.savetxt("best_kin_specs.txt", new_kin_data_out, delimiter = "\t", fmt = "%1.8f")
        with open("kin_specs_history.txt", 'a') as f:
            for item in new_kin_data:
                f.write("%1.8f\t" % item)
            f.write("\n")
        # Write the variance to "smallest var .txt"
        mini_data_out = np.full(1,mini)
        np.savetxt("smallest_var.txt", mini_data_out, delimiter = "\t", fmt = "%1.8f")
        # Save the current generation associated with the new minimum
        mini_data_out[0] = int(sys.argv[2])
        np.savetxt("gen_of_smallest_var.txt", mini_data_out, delimiter = "\t", fmt = "%d")
        print( "mini_difference pct = " + str((previous_mini-mini)/previous_mini));
        if (previous_mini-mini)/previous_mini > noise_range:
            count_unsuccessfull = 0
    #        np.savetxt("no_reduction_count.txt",0)
        else:
            count_unsuccessfull += 1
    #        np.savetxt("no_reduction_count.txt",count_unsuccessfull)
    else:
        count_unsuccessfull +=1

    count_unsuccessfull_out[0] = count_unsuccessfull    
    np.savetxt("no_reduction_count.txt",count_unsuccessfull_out)

    mini_posi_and_var = np.full(1,mini_nr)
    np.savetxt("best_candidate.txt", mini_posi_and_var, delimiter = "\t", fmt = "%1i")




