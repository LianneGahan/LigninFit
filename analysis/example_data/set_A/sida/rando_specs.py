from numpy import full,savetxt
import sys
import numpy as np
import sys
import random
import os


def rnum():#return random number between -1 and 1
    val = random.gauss(0,0.25)
    while val < -1 or val > 1:
        val = random.gauss(0,0.25)
    return val

#argv[1] = generation
#argv[2] = N_runs
#argv[3] = delta
#argv[4] = params_to_change
#argv[5] = parentfolder


def gradient():
    print("CALCULATING GRADIENT")

    #Maximum values of each parameter
    max_vals_kin = np.loadtxt("max_kin_vals.txt")

    #Relative minimum values of each parameter
    rel_min_vals_kin = np.loadtxt("min_kin_vals.txt")
    rel_min_vals_kin[:] /= max_vals_kin[:]

    #Generation N-2
    old_var = np.loadtxt("old_smallest_var.txt");
    old_kin_data = np.loadtxt("old_best_kin_specs.txt")#//
    old_kin_data[:] /= max_vals_kin[:]

    #Generation N-1
    new_var = np.loadtxt("smallest_var.txt");
    new_kin_data = np.loadtxt("best_kin_specs.txt")
    new_kin_data[:] /= max_vals_kin[:]

    print("Relative mins: ", rel_min_vals_kin)
    print("old kins: ", old_kin_data)
    print("new kins: ", new_kin_data)

    #Denotes which parameters should be varied
    variation_bools_kin = np.loadtxt("kin_params_to_randomize.txt", dtype = int)#Contains boolean values for which parameters are varied

    N_kin_vals = new_kin_data.size


    kin_gradient = [0. for i in range(N_kin_vals)]
    for i in range(N_kin_vals):
        #print(new_var, old_var, new_kin_data[i], old_kin_data[i])
        if new_kin_data[i] - old_kin_data[i] != 0:
            kin_gradient[i] = (new_var - old_var) / (new_kin_data[i] - old_kin_data[i])
            #print("Gradient value ", i, kin_gradient[i])

    print("Gradient values ", kin_gradient)
    out_kin_data = np.full(N_kin_vals,0.)


    for j in range(int(N_runs)):
        step_size = delta;
        for i in range(N_kin_vals):
            if variation_bools_kin[i] == 1:
                value_to_set = new_kin_data[i]*kin_gradient[i]
                print(kin_gradient[i])
                #print("Delta is: ", delta)
                #print("Value to set is: ", value_to_set, " for variable ", i)
                #print("relative minimum is ", rel_min_vals_kin[i])
                if value_to_set > rel_min_vals_kin[i] and value_to_set <= 1.:
                    out_kin_data[i] = value_to_set
                    print("Changing kin param %d from %1.4f to %1.4f" % (i, new_kin_data[i], value_to_set))

                else:
                    print("Not changing kin param %d from %1.4f to %1.4f" % (i, new_kin_data[i], value_to_set))
                    out_kin_data[i] = new_kin_data[i]
            else:
                out_kin_data[i] = new_kin_data[i]
        with open("rando_specs_history.txt", 'a') as f:
            f.write("%1.8f\t" % int(generation))
            for item in out_kin_data:
                f.write("%1.8f\t" % item)
            f.write("\n")
        with open("family_" + parent_folder + "/Generation_" + str(generation) + "/Run_" + str(j+1) + "/Params/kinetic_parameters.txt","w") as f:
            for k, data_point in enumerate(out_kin_data):
                val = data_point * max_vals_kin[k]
                f.write("%f\t" % val)



################################################################################################################


def random_sampling():

    #Maximum values of each parameter
    max_vals_kin = np.loadtxt("max_kin_vals.txt")

    #Relative minimum values of each parameter
    rel_min_vals_kin = np.loadtxt("min_kin_vals.txt")
    rel_min_vals_kin[:] /= max_vals_kin[:] 

    #Generation N-1
    new_var = np.loadtxt("smallest_var.txt");
    new_kin_data = np.loadtxt("best_kin_specs.txt") / max_vals_kin[:]

    #Denotes which parameters should be varied
    variation_bools_kin = np.loadtxt("kin_params_to_randomize.txt", dtype = int)#Contains boolean values for which parameters are varied

    N_kin_vals = new_kin_data.size

    out_kin_data = np.full(N_kin_vals,0.)

    for j in range(int(N_runs)):
        for i in range(N_kin_vals):
            step_size = rnum()*delta;
            if variation_bools_kin[i] == 1:
                value_to_set = new_kin_data[i] + step_size;
                if value_to_set > rel_min_vals_kin[i] and value_to_set <= 1. and value_to_set >= 0.:
                    print("Changing kin param %d from %1.4f to %1.4f" % (i,new_kin_data[i],value_to_set))
                    out_kin_data[i] = value_to_set
                else:
                    print("Not changing kin param %d from %1.4f to %1.4f" % (i,new_kin_data[i],value_to_set))
                    out_kin_data[i] = new_kin_data[i]
            else:
                out_kin_data[i] = new_kin_data[i]
        with open("family_" + parent_folder + "/Generation_" + str(generation) + "/Run_" + str(j+1) + "/Params/kinetic_parameters.txt","w") as f:
            for k,data_point in enumerate(out_kin_data):
                val = data_point * max_vals_kin[k]
                #print("VALUE ", k, val)
                f.write("%f\t" % val)
        with open("rando_specs_history.txt", 'a') as f:
            f.write("%1.8f\t" % int(generation))
            for item in out_kin_data:
                f.write("%1.8f\t" % item)
            f.write("\n")


        
################################################################################################################


generation=sys.argv[1]
N_runs=int(sys.argv[2])
delta=float(sys.argv[3])
#current_params =sys.argv[4]
parent_folder = sys.argv[4]

N_unsuccessful_tries = np.loadtxt("no_reduction_count.txt")
if N_unsuccessful_tries > 2:
    delta += delta * N_unsuccessful_tries/10.


print("Delta = " + str(delta))


print("In rando_specs.py: N_unsuccessful_tries = " + str(N_unsuccessful_tries))

print("Generation = " + str(generation))
if N_unsuccessful_tries < 1 and int(generation) != 1:   
    print("Gradient param update");
    gradient()


else:

    print("Random param update (no gradient used)"); ##We also need a term here to use the original conditions set by the grid search
    random_sampling()





