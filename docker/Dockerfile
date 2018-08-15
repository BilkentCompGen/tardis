FROM ubuntu:18.04

RUN apt-get update -y --fix-missing
RUN apt-get install git make gcc zlib1g-dev libbz2-dev liblzma-dev -y

RUN mkdir /tardis
RUN mkdir /input
RUN mkdir /output
WORKDIR /tardis

RUN git clone --recursive https://github.com/BilkentCompGen/tardis.git /tardis
RUN make libs && make

ENV PATH="/tardis/bin:${PATH}"
VOLUME /input
VOLUME /output
ENTRYPOINT ["/tardis/tardis"]
