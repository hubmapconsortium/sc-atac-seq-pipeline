FROM ubuntu:18.04

RUN apt-get update
RUN apt-get install -y python3-pip python-pip samtools tabix
RUN pip3 install sinto

COPY sort_index_frag.py /opt/
