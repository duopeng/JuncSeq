# JuncSeq
Extract splice junction sequences (from Ensembl format genomes)  
  
![image](https://user-images.githubusercontent.com/4129442/220262088-47984630-fc4c-4448-914d-69da66d0834e.png)



Make Weblogo plots  
  
![image](https://user-images.githubusercontent.com/4129442/220261597-416a3bf5-1acd-49c9-89dc-416be11e791f.png)
![image](https://user-images.githubusercontent.com/4129442/220261672-3472b0ea-36fb-4311-8a31-71744907963d.png)  

![image](https://user-images.githubusercontent.com/4129442/220261624-1336c745-a957-4178-af1f-93996ef9dace.png)
![image](https://user-images.githubusercontent.com/4129442/220261686-5793cb04-193e-4fad-bbce-2f2ba8b0e3d8.png)  




# Usage

### Clone the repository
```
git clone https://github.com/duopeng/JuncSeq
```
### Go the repository directory
```
cd JuncSeq
```
### Create conda environment and activate it
```
conda create -y -n JuncSeq python=3.10
conda activate JuncSeq
```
### Install required packages
```
pip install -r requirements.txt
```
## Run tests to verify installation
```
python run_tests.py
```
You should see the following output:  
  
![image](https://user-images.githubusercontent.com/4129442/220225526-c46da6d6-9a04-4d1f-95b2-0b8615e13a5b.png)

## Use the following notebook to extract splice junction sequences and make Weblogo plots
```
JuncSeq.ipynb
```
Please make sure the "GenomeBuild" and "species" match, see the mouse example below:  

![image](https://user-images.githubusercontent.com/4129442/220231238-069235ad-6561-4fe5-a917-39e98173fb4e.png)



