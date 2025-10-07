##############################################
# FUNCTION: sacct_summary
# DESCRIPTION:
#   Retrieve and summarize SLURM accounting information for a job.
#   Extracts key job metrics (runtime, CPU usage, memory, nodes, etc.)
#   and outputs a clean summary in parsable tabular format.
##############################################
sacct_summary() {
    if [ -z "$1" ]; then
        echo "Usage: sacct_summary <JobID>"
        return 1
    fi

    job_id="$1"

    sacct -j "$job_id" --format=JobID,JobName,State,Elapsed,TotalCPU,ReqCPUS,Timelimit,MaxRSS,NNodes,NTasks,ReqMem --parsable2 | awk -F '|' '
    NR==1 { print "JobID|JobName|State|Elapsed|Timelimit|TotalCPU|ReqCPUS|UsedCPUS|MaxRSS|NNodes|NTasks|ReqMem"; next }
    $2 !~ /batch|extern/ { jobname=$2; reqmem=$11; timelimit=$7; reqcpus=$6 }  # Store main job info
    $2 == "batch" {
        sub(/\.batch/, "", $1)

        mem_kb = ($8+0)
        mem_gb = (mem_kb / (1024^2))
        mem_rec = (mem_gb * 1.2) * reqcpus

        split($4, elapsed, ":")
        if (length(elapsed) == 3) { elapsed_sec = elapsed[1] * 3600 + elapsed[2] * 60 + elapsed[3] }
        else if (length(elapsed) == 2) { elapsed_sec = elapsed[1] * 60 + elapsed[2] }
        else { elapsed_sec = elapsed[1] }

        split($5, totalcpu, ":")
        if (length(totalcpu) == 3) { totalcpu_sec = totalcpu[1] * 3600 + totalcpu[2] * 60 + totalcpu[3] }
        else if (length(totalcpu) == 2) { totalcpu_sec = totalcpu[1] * 60 + totalcpu[2] }
        else { totalcpu_sec = totalcpu[1] }

        if (elapsed_sec > 0) {
            used_cpus = totalcpu_sec / elapsed_sec
        } else {
            used_cpus = 0
        }

        split(timelimit, reqtime_arr, ":")
        if (length(reqtime_arr) == 3) { req_sec = reqtime_arr[1] * 3600 + reqtime_arr[2] * 60 + reqtime_arr[3] }
        else if (length(reqtime_arr) == 2) { req_sec = reqtime_arr[1] * 60 + reqtime_arr[2] }
        else { req_sec = reqtime_arr[1] }

        rec_sec = int(elapsed_sec * 1.2)
        rec_hh = int(rec_sec / 3600)
        rec_mm = int((rec_sec % 3600) / 60)
        rec_ss = rec_sec % 60
        recommended_time = sprintf("%02d:%02d:%02d", rec_hh, rec_mm, rec_ss)

        printf "%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n", $1, jobname, $3, $4, timelimit, $5, reqcpus, $8, $9, $10, reqmem
    }'
}