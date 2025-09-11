# Use a Miniconda base image
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy environment.yml into the container
COPY environment.yml .

# Create the conda environment
RUN conda env create -f environment.yml

# Make sure conda is initialized
SHELL ["conda", "run", "-n", "DGM-env", "/bin/bash", "-c"]

# (Optional) set default environment
ENV PATH /opt/conda/envs/DGM-env/bin:$PATH

# Copy your project files
COPY . . 


