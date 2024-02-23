FROM public.ecr.aws/ubuntu/ubuntu:20.04_stable

# install AWS CLI and Python 3.
# - aws CLI is required by nextflow to copy the input and output data.
# - salmon and fastqc are used in the workflow .

ENV FASTQC_URL https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
ENV SALMON_URL https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y \
    default-jre \
    python3-pip \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# install awscli
RUN pip3 install awscli --upgrade
RUN ln -s /usr/bin/python3 /usr/local/bin/python

# install FastQC
RUN wget -O /tmp/z.$$ ${FASTQC_URL} && (cd /usr/local/ && unzip /tmp/z.$$) && rm /tmp/z.$$
RUN ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc

# install salmon
RUN wget -O /tmp/salmon.tar.gz ${SALMON_URL} && (cd /usr/local/ && tar xzf /tmp/salmon.tar.gz) && rm /tmp/salmon.tar.gz
RUN ln -s /usr/local/salmon-latest_linux_x86_64/bin/salmon /usr/local/bin/salmon

