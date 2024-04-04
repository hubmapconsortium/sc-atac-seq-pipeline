FROM pypy:3

## Python packages
COPY docker/barcode-adj/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

COPY bin /opt
COPY data/* /opt/

WORKDIR /opt
