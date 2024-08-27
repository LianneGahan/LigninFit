import numpy as np

N_data_points = np.loadtxt("latest/Params/kinetic_parameters.txt").size;
out_file = np.full((1,N_data_points),0.);
np.savetxt("current_gradient.txt", out_file, delimiter = "\t", fmt = "%1.8f")


N_data_points = np.loadtxt("latest/Params/kinetic_parameters.txt").size;
out_file = np.full((1,N_data_points),0.);
np.savetxt("old_best_kin_specs.txt", out_file, delimiter = "\t", fmt = "%1.8f")

N_data_points = np.loadtxt("old_smallest_var.txt").size;
out_file = np.full((1,N_data_points),12.);
np.savetxt("old_smallest_var.txt", out_file, delimiter = "\t", fmt = "%1.8f")

#N_data_points = np.loadtxt("smallest_var.txt").size;
#out_file = np.full((1,N_data_points),1.);
#np.savetxt("smallest_var.txt", out_file, delimiter = "\t", fmt = "%1.8f")
