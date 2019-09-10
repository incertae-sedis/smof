FROM alpine:3.8

RUN apk update
RUN apk add --no-cache \
    emacs \
    git \
    python3

#RUN git clone https://github.com/incertae-sedis/smof.git smof-git
#RUN ln -s /smof-git/smof.py /smof
RUN pip3 install smof
RUN mkdir inout
WORKDIR /inout
ENV PATH=/:/smof-git/:$PATH
CMD /smof ; \
    echo -e "\nExample Docker Run: docker run -v \${PWD}:/inout/ smof:latest /smof grep Iowa fasta.fa > Iowa.fa"

LABEL author="Jennifer Chang"
LABEL last-update="2019-09-10"
