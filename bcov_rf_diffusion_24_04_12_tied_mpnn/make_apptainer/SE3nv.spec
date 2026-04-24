Bootstrap: docker
From: ubuntu
IncludeCmd: yes

%setup
rsync -a --no-g --no-o ${LOGOS_PATH%/}/bcov_rf_diffusion_24_04_12_tied_mpnn/SE3Transformer/ $APPTAINER_ROOTFS/SE3Transformer/

%files
/etc/localtime
/etc/hosts
/etc/apt/sources.list
/storage/d1/users/ruizhi/softwares/Miniconda3-latest-Linux-x86_64.sh /opt/miniconda.sh

%post
# Switch shell to bash
rm /bin/sh; ln -s /bin/bash /bin/sh

# Common symlinks
ln -s /net/databases /databases
ln -s /net/software /software
ln -s /home /mnt/home
ln -s /projects /mnt/projects
ln -s /net /mnt/net

apt-get update
# required X libs
apt-get install -y libx11-6 libxau6 libxext6 libxrender1
#  git
apt-get install -y git
apt-get clean

# Install conda
bash /opt/miniconda.sh -b -u -p /usr

# Install conda/pip packages
# conda update conda

# Initialize conda for bash
eval "$(conda shell.bash hook)"

conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create environment with python 3.9
CONDA_OVERRIDE_CUDA='12.4' conda create  -n rfdiffusion python=3.9 -y

# Activate the environment
conda activate rfdiffusion

# open-babel
conda install  -n rfdiffusion \
   openbabel \
   -c conda-forge \
   --solver=classic \
   -y

# pytorch + dependancies
conda install  -n rfdiffusion \
   -c nvidia \
   -c pytorch \
   -c pyg \
   -c dglteam/label/cu124 \
   -c https://conda.rosettacommons.org \
   pip \
   ipython \
   ipykernel \
   numpy \
   pandas \
   seaborn \
   matplotlib \
   jupyterlab \
   pytorch \
   pytorch-cuda=12.4 \
   pyg \
   pyrosetta \
   -y

pip install --no-cache-dir --no-dependencies dgl -f https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html
# pip extras
pip install e3nn \
   omegaconf \
   hydra-core \
   pyrsistent \
   opt_einsum \
   sympy \
   icecream \
   wandb \
   deepdiff \
   scikit-learn \
   assertpy

# SE3 transformer
pip install /SE3Transformer/

# Clean up
conda clean -a -y
apt-get -y purge build-essential git wget
apt-get -y autoremove
apt-get clean
rm /opt/miniconda.sh

%environment
export PATH=$PATH:/usr/local/cuda/bin
# Activate conda environment
source /usr/etc/profile.d/conda.sh
conda activate rfdiffusion

%runscript
source /usr/etc/profile.d/conda.sh
conda activate rfdiffusion
exec python "$@"

%help
SE3nv environment for running RF-diffusion, etc.
