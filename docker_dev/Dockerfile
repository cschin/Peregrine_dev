FROM continuumio/miniconda3
RUN apt-get update
RUN apt-get install -y build-essential zlib1g zlib1g-dev
RUN apt-get install -y libbz2-dev
RUN mkdir /opt/build
COPY src.tgz /opt/build
COPY install_with_conda.sh /opt/build
RUN cd /opt/build; tar zxvf src.tgz; bash install_with_conda.sh
RUN . /opt/conda/bin/activate; conda clean --all
RUN apt-get update
RUN apt-get install -y parallel time
RUN . /opt/conda/bin/activate; conda activate peregrine;
RUN apt-get install -y make
RUN mkdir /opt/licenses
COPY LICENSE /opt/licenses/LICENSE
COPY LICENSE.falcon /opt/licenses/LICENSE.falcon
COPY LICENSE.minimap2 /opt/licenses/LICENSE.minimap2
RUN mkdir /opt/test
COPY test/Makefile /opt/test
COPY test/run_test.sh /opt/test
COPY test/simulate_reads.py /opt/test
COPY bashrc /root/.bashrc
COPY entry.sh /opt/
COPY entry_dev.sh /opt/
COPY jupyterlab.sh /opt/
WORKDIR /opt/test
ENTRYPOINT ["/opt/entry_dev.sh"]
