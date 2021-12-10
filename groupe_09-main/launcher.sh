### launcher
### instructions required to prepare an IFB VM (16 CPU, 64Go)

<<COMMENT #To allow parallelisation on file to download :
We wanted to parallelise the downloading rule to gain time.
That means the iteration is the downloading rule itself (i.e. the rule is generalised) but that implies NO FOR LOOP (on the sample) inside the rule.
Thus the rule needs as input the names of sample (as {sample}). Alas their names are inside a list, not actual file. Thus the need to create empty files.
This is done with the touch command.
COMMENT

listeDownload=('SRR628585' 'SRR628586' 'SRR628587' 'SRR628588' 'SRR628582' 'SRR628584' 'SRR628583' 'SRR628589' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'MT')
for i in ${listeDownload[*]}; do
        touch ${i} #create empty files
done

# Install R
sudo add-apt-repository "deb https://stat.ethz.ch/CRAN/bin/linux/ubuntu bionic-cran35/"
sudo apt install r-base r-base-dev -y
sudo apt-get install libcurl4-openssl-dev -y
sudo apt-get install libxml2-dev -y

# launching pipeline analysis
snakemake -s pipeline.smk --cores 2 --use-singularity