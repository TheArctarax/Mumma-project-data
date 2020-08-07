import numpy as np

job_name = "bilby_run_distance"
var_name = "distance"
exp_name = "lambda_bimodality"
start_val = 110
end_val = 490
val_spacing = 10

job_array = np.linspace(start_val, end_val, int((end_val - start_val) / val_spacing))
number_of_jobs = len(job_array)

with open("bilby_runs_test.dag", "w", encoding="utf-8") as f:
    for job_number in range(1, number_of_jobs + 1):
        f.write("JOB {}_{} bilby_runs.sub\n".format(job_name, str(job_number)))
        f.write("VARS {}_{}".format(job_name, str(job_number)))
        f.write(" {}=".format(var_name.upper()))
        f.write(
            "'{}' OUTDIR='/home/darin.mumma/public_html".format(
                str(round(job_array[job_number - 1]))
            )
        )
        f.write("/bilby_noiseless/{}'".format(exp_name))
        f.write(" LABEL='{}_".format(var_name.lower()))
        f.write("{}'".format(str(round(job_array[job_number - 1]))))
        f.write(" EXPERIMENT='{}'".format(exp_name))
        f.write(" RUN='{}'\n\n".format(str(job_number)))
