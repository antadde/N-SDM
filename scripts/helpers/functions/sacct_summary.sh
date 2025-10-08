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

    sacct -j "$job_id" \
        --format=JobID,JobName,State,Elapsed,TotalCPU,ReqCPUS,Timelimit,MaxRSS,NNodes,NTasks,ReqMem \
        --parsable2 | \
    awk -F '|' '
    NR==1 {
        print "JobID|JobName|State|Elapsed|Timelimit|ReqCPUS|UsedCPUS|NNodes|NTasks|MaxRSS|ReqMem"
        next
    }

    # Skip extern/noise jobs
    $2 !~ /batch|extern/ {
        jobname = $2
        reqmem = $11
        timelimit = $7
        reqcpus = $6
    }

    # Process main batch lines
    $2 == "batch" && $4 ~ /:/ {
        sub(/\.batch/, "", $1)

        # --- Convert Elapsed and TotalCPU to seconds ---
        split($4, elapsed, ":")
        elapsed_sec = (length(elapsed) == 3) ? elapsed[1]*3600 + elapsed[2]*60 + elapsed[3] :
                      (length(elapsed) == 2) ? elapsed[1]*60 + elapsed[2] : elapsed[1]

        split($5, totalcpu, ":")
        totalcpu_sec = (length(totalcpu) == 3) ? totalcpu[1]*3600 + totalcpu[2]*60 + totalcpu[3] :
                       (length(totalcpu) == 2) ? totalcpu[1]*60 + totalcpu[2] : totalcpu[1]

        used_cpus = (elapsed_sec > 0) ? totalcpu_sec / elapsed_sec : 0

        # --- Convert MaxRSS to GB ---
        mem = $8
        mem_val = mem
        if (mem ~ /K$/) { sub(/K$/, "", mem_val); mem_gb = mem_val / (1024^2) }
        else if (mem ~ /M$/) { sub(/M$/, "", mem_val); mem_gb = mem_val / 1024 }
        else if (mem ~ /G$/) { sub(/G$/, "", mem_val); mem_gb = mem_val + 0 }
        else { mem_gb = 0 }
        mem_gb = sprintf("%.1fG", mem_gb)

        printf "%s|%s|%s|%s|%s|%s|%.2f|%s|%s|%s|%s\n",
            $1, jobname, $3, $4, timelimit, reqcpus, used_cpus, $9, $10, mem_gb, reqmem
    }'
}