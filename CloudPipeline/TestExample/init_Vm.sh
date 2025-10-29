########### instal docker
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "${UBUNTU_CODENAME:-$VERSION_CODENAME}") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
# To install the latest version, run:
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

#Create the docker group.
sudo groupadd docker
#Add your user to the docker group.
sudo usermod -aG docker $USER
# activate the changes to groups
newgrp docker
# Verify that you can run docker commands without sudo
docker run hello-world

################# install Nextflow 
install zip
sudo apt install zip
#Install SDKMAN:
curl -s https://get.sdkman.io | bash
## !!! open new terminal!!!
#Install Java:
sdk install java 17.0.10-tem
java -version
# Install Nextflow
curl -s https://get.nextflow.io | bash
# make it runnable
chmod +x nextflow
mkdir bin
mv nextflow bin
echo "export PATH=${HOME}/bin:\$PATH" >> ${HOME}/.bashrc
# !! open new terminal!!"
nextflow info

# ################# install aws cli
# ### install aws cli
# cd bin
# curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
# unzip awscliv2.zip
# sudo ./aws/install

# ################# install gsutil cli
# curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-linux-x86_64.tar.gz
# tar -xf google-cloud-cli-linux-x86_64.tar.gz
# ./google-cloud-sdk/install.sh
# ./google-cloud-sdk/bin/gcloud init