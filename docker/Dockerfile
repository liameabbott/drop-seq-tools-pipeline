FROM openjdk:8-slim

ARG drop_seq_tools_version=2.4.0
ARG picard_version=2.20.5
ARG star_version=2.7.5c
ARG samtools_version=1.11

WORKDIR /working

RUN apt-get update \
	&& apt-get install -y \
		wget \
		unzip \
		gcc \
		bzip2 \
		libncurses5-dev \
		zlib1g-dev \
		libbz2-dev \
		liblzma-dev \
		make \
		file \
	&& rm -rf /var/lib/apt/lists/* \
	&& mkdir -p /working/software/bin/ \
	&& cd /working/software/ \
	# install Drop-seq tools
	&& wget https://github.com/broadinstitute/Drop-seq/releases/download/v${drop_seq_tools_version}/Drop-seq_tools-${drop_seq_tools_version}.zip \
	&& unzip Drop-seq_tools-${drop_seq_tools_version}.zip  \
	&& mv Drop-seq_tools-${drop_seq_tools_version}/ /working/software/drop-seq-tools \
	&& rm -r Drop-seq_tools-${drop_seq_tools_version}.zip \
	&& ln -s /working/software/drop-seq-tools/* /working/software/bin/ \
	&& ln -s /working/software/drop-seq-tools/jar/lib/picard-${picard_version}.jar /working/software/picard.jar \
	# install AWS CLI
	&& wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip \
	&& unzip awscli-exe-linux-x86_64.zip \
	&& ./aws/install -i /working/software/aws-cli \
	&& rm -r awscli-exe-linux-x86_64.zip aws \
	&& ln -s /working/software/aws-cli/v2/current/bin/* /working/software/bin/ \
	# install STAR
	&& wget https://github.com/alexdobin/STAR/archive/${star_version}.zip \
	&& unzip ${star_version}.zip \
	&& mv /working/software/STAR-${star_version}/ /working/software/star/ \
	&& rm -r ${star_version}.zip \
	&& ln -s /working/software/star/bin/Linux_x86_64/STAR /working/software/bin/ \
	# install samtools
	&& wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
	&& tar -xjvf samtools-${samtools_version}.tar.bz2 \
	&& cd samtools-${samtools_version} \
	&& ./configure --prefix=/working/software/samtools \
	&& make \
	&& make install \
	&& cd /working/software \
	&& rm -r samtools-${samtools_version}.tar.bz2 samtools-${samtools_version} \
	&& ln -s /working/software/samtools/bin/* /working/software/bin/

ENV PATH="/working/software/bin:$PATH"
