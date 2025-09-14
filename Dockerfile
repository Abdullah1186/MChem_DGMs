# Use a Miniconda base image
FROM continuumio/miniconda3

# Set working directory
WORKDIR /MChem_DGMs

# Copy environment.yml into the container
COPY environment.yml .

# Create the conda environment
RUN conda env create -f environment.yml

# Make sure conda is initialized
SHELL ["conda", "run", "-n", "DGM-env", "/bin/bash", "-c"]

# (Optional) set default environment
ENV PATH=/opt/conda/envs/DGM-env/bin:$PATH

# Copy your project files
COPY . . 

# Make a directory for your data
RUN mkdir -p /MChem_DGMs/analysis/Databases \
    && wget https://moleculardatabases.s3.eu-west-2.amazonaws.com/databases.tar.gz -O /tmp/databases.tar.gz \
    && tar -xvzf /tmp/databases.tar.gz -C /MChem_DGMs/analysis/Databases \
    && rm /tmp/databases.tar.gz



