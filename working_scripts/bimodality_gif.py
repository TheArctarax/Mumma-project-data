import imageio
import numpy as np

images = []
filenames = []

start_val=110
end_val=490
val_spacing=10

job_array = np.arange(start_val, end_val + val_spacing, val_spacing)
number_of_jobs = len(job_array)

for job_number in range(1, number_of_jobs+1):
    filenames.append("./processed_{}.png".format(job_array[job_number-1]))

for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('./lambda_bimodality.gif', images, fps=2)

