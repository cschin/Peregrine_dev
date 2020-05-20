FROM continuumio/miniconda3
RUN apt-get update
RUN apt-get install -y build-essential zlib1g zlib1g-dev
RUN apt-get install -y libbz2-dev
RUN mkdir /opt/build
COPY src.tgz /opt/build
COPY install_with_conda.sh /opt/build
RUN cd /opt/build; tar zxvf src.tgz; bash install_with_conda.sh
RUN . /opt/conda/bin/activate; conda clean -all

FROM continuumio/miniconda3
COPY --from=0 /opt/conda /opt/conda
RUN apt-get update
RUN apt-get install -y parallel time
RUN . /opt/conda/bin/activate; conda activate peregrine;
RUN apt-get install -y make
RUN mkdir /opt/licenses
COPY LICENSE /opt/licenses/LICENSE
COPY LICENSE.falcon /opt/licenses/LICENSE.falcon
COPY LICENSE.minimap2 /opt/licenses/LICENSE.minimap2
RUN mkdir -p /wd
COPY bashrc /root/.bashrc
WORKDIR /wd
RUN cd /opt/ && wget -q https://github.com/lh3/dipcall/releases/download/v0.1/dipcall-0.1_x64-linux.tar.bz2
RUN cd /opt/ && tar -jxf dipcall-0.1_x64-linux.tar.bz2 --no-same-owner
COPY generate_cns.py /opt/generate_cns.py
COPY entry.sh /opt/entry.sh 

ENTRYPOINT ["/opt/entry.sh"]
