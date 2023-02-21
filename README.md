# JuncSeq
Extract splice junction sequences (from Ensembl format genomes)  
  
![image](https://user-images.githubusercontent.com/4129442/220227109-fd0b8f78-338a-4c37-bc71-b63e24f5bd6c.png)


Make Weblogo plots  
  
![image](https://user-images.githubusercontent.com/4129442/220224109-0fbda8bb-485c-4ce4-a080-a918380298a4.png)
![image](https://user-images.githubusercontent.com/4129442/220224160-ced5e3c6-4c59-48f0-b839-c3b3f3ccfd32.png)


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


