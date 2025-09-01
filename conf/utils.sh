#!/bin/bash

#############################################################################
# Log functions
#############################################################################
# Function for logging message at specific time
log_message() {
    echo "[$(date +"%F %T")] $1"
}

export -f log_message

log_error() {
    local timestamp=$(date +'%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ERROR: $*" >&2 | tee -a "${LOG_FILE:-/dev/stderr}"
}


log_warning() {
    local timestamp=$(date +'%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] WARNING: $*" | tee -a "${LOG_FILE:-/dev/stdout}"
}

# Function to calculate and display elapsed time
log_time_usage() {
    local start_time=$1
    local end_time=${2:-$(date +%s)}  # Use current time if end_time not provided
    
    local total_time=$((end_time - start_time))
    local hours=$((total_time / 3600))
    local minutes=$(( (total_time % 3600) / 60 ))
    local seconds=$((total_time % 60))
    
    echo "[$(date +"%F %T")] Total runtime: ${hours}h:${minutes}m:${seconds}s"
}

export -f log_time_usage

# Function to check if a command executed successfully
check_command() {
    if [ $? -ne 0 ]; then
        log_message "ERROR: $1 failed"
        return 1
    fi
}

export -f check_command

# Function to check if a file exists
validate_inputs() {
    local input_file=$1
    [[ ! -f "$input_file" ]] && { 
        log_message "ERROR: Input file not found: $input_file"
        return 1
    }
}