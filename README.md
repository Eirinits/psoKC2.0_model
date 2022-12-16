# The extended Psoriatic Keratinocyte model (PsoKC 2.0)
A repository for the extended logical model of psoriatic keratinocytes, based on the [PsoKC model by Tsirvouli et al.](https://www.sciencedirect.com/science/article/pii/S258900422101422X)

The model was built by Eir Aker and refined by Eirini Tsirvouli. Both contributed to its analysis. 

This repository contains ipynb files and psoriatic keratinocyte model-related files as described in paper "Patient-specific logical models replicate phenotype responses to psoriatic and anti-psoriatic stimuli". 

## How to run this notebook?

### With docker

1. Install the [colomoto docker image](https://github.com/colomoto/colomoto-docker). 
The model presented in the paper was tested with the version ```2020-01-24```.
2. Install [jupyter notebook](http://jupyter.org/).
3. Clone the repository: ```https://github.com/Eirinits/psoKC2.0_model.git```. 
4. Go to the project repository and launch the notebook: ```colomoto-docker --bind .```
5. The notebook will be available in your web browser at ```http://localhost:8888/```
