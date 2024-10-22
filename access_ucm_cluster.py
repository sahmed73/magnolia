# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Sep 12 01:06:01 2024
"""

import paramiko
import tempfile

def connect():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    
    # Connect to the HPC via SSH
    hostname = 'login.rc.ucmerced.edu'
    username = 'sahmed73'
    password = '$i1oveMyself$'  # Secure password handling recommended

    try:
        ssh.connect(hostname, username=username, password=password)
        print('Connected to the UCM HPC cloud!!')
        return ssh
    except paramiko.AuthenticationException:
        print("Authentication failed, please verify your credentials.")
        return None
    except paramiko.SSHException as e:
        print(f"Unable to establish SSH connection: {e}")
        return None
    except Exception as e:
        print(f"Exception in connecting to the server: {e}")
        return None

def local_copy(remote_file):
    ssh = connect()
    if ssh is None:
        print("Connection failed, cannot proceed.")
        return None

    try:
        # Open SFTP session to transfer the file
        sftp = ssh.open_sftp()

        # Create a temporary file
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        local_temp_path = temp_file.name
        temp_file.close()

        # Download the file to the temporary location
        # print(f"Downloading {remote_file} to {local_temp_path}...")
        sftp.get(remote_file, local_temp_path)
        print(f"File downloaded to {local_temp_path}!")

        # Close the SFTP session and SSH connection
        sftp.close()
        ssh.close()

        # Return the local temporary file path for further use
        return local_temp_path

    except Exception as e:
        print(f"Error during file transfer: {e}")
        if ssh:
            ssh.close()
        return None
