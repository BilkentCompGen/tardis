FROM ubuntu:18.04

RUN apt-get update -y --fix-missing
RUN apt-get upgrade -y

RUN apt-get install git make gcc zlib1g-dev libbz2-dev liblzma-dev gnuplot-nox -y

RUN mkdir /tardis
RUN mkdir /input
RUN mkdir /output
WORKDIR /tardis

RUN git clone --recursive https://github.com/BilkentCompGen/tardis.git /tardis
RUN make libs && make

RUN mkdir /mrfast
WORKDIR /mrfast
RUN git clone https://github.com/BilkentCompGen/mrfast.git /mrfast
RUN make

WORKDIR /tardis

ENV PATH="/tardis:/mrfast:${PATH}"
RUN apt-get remove git -y
RUN apt-get autoremove -y
VOLUME /input
VOLUME /output
ENTRYPOINT ["/tardis/tardis"]
