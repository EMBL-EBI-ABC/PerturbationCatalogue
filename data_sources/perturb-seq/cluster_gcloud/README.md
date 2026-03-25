# Google Cloud SDK Installation on Slurm Cluster (Non-Root)

This guide explains how to install and authenticate the `gcloud` and `gsutil` commands on a Slurm cluster where you do not have root privileges.

All commands must be executed inside an interactive session.

## 1. Download and Install Google Cloud SDK

```bash
# Create the software directory if it doesn't exist
mkdir -p "$HPS_PATH/software"

# Download the Linux 64-bit archive
curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-468.0.0-linux-x86_64.tar.gz

# Extract the archive into $HPS_PATH/software
tar -xvf google-cloud-sdk-468.0.0-linux-x86_64.tar.gz -C "$HPS_PATH/software"

# Run the installation script
# This script will ask if you want to update your PATH in .bashrc (say yes)
"$HPS_PATH/software/google-cloud-sdk/install.sh"
```

## 2. Update Environment

After installation, you need to refresh your shell to include the `gcloud` binaries in your `PATH`.

```bash
source ~/.bashrc
```

## 3. Authenticate

On a remote cluster, you must use the `--no-launch-browser` flag to authenticate, as there is no graphical browser available.

```bash
gcloud auth login --no-launch-browser
```

1.  Copy the URL provided in the terminal into your local machine's browser.
2.  Log in with your Google account.
3.  Copy the authorization code provided.
4.  Paste it back into the terminal on the cluster.

## 4. Verify Installation

Check if the commands are working:

```bash
gcloud version
gsutil version
```

## 5. (Optional) Cleanup

You can remove the downloaded archive after installation:

```bash
rm google-cloud-sdk-468.0.0-linux-x86_64.tar.gz
```
