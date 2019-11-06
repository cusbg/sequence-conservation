# Sequence conservation calculator

## Quick start with Docker
The following lines demonstrates how to setup Sequence conservation on 
local computer using Docker. 

### Installation

Clone this repository:
```
git clone https://github.com/skodapetr/sequence-conservation.git
```

Enter into the new folder and run Docker build:
```
cd sequence-conservation
docker build \
    -t conservation \
    --build-arg USER=$(id -u ${USER}) \
    .
```

Run Docker image with mapping data directory into ```/data/conservation``` 
on host: 
```
docker run \
    -v /data/conservation:/data/conservation \
    -u $(id -u ${USER}):$(id -g ${USER}) \
    -it conservation
```

In the next step we need to prepare a database for BLAST. There
are two options, the database can be created locally on downloaded.
Use one of the following command to prepare the database:
```
./download_database.sh
```
or
```
./prepare_database.sh
```

### Demonstration

We can showcase use of the script by calculating conservation for ```2SRC```.
First we enter the ```/data/conservation``` directory and download the FASTA file 
from PDB, next we run the Python script that computes the conservation. 
```
cd /data/conservation
curl "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=2SRC&compressionType=uncompressed" -o 2SRC.fasta
python3 /opt/conservation/calculate_conservation.py --input ./2SRC.fasta --output ./2SRC.json
```

You can close the Docker container by typing ```exit```.
When the Docker container is terminated all data outside volumes
are lost. As we mapped the ```/data/conservation``` directory, the data 
in this directory are preserved and available from the host computer.

Another advantage of mapping of ```/data/conservation``` is 
that the BLAST database is preserved and does not need to be 
rebuild each time a new Docker container is started.

## Credits
The initial content of this repository was heavily inspired by [calc_protein_conservation](https://github.com/jendelel/calc_protein_conservation).
