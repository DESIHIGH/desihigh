FROM ubuntu:latest

ENV TZ=US

RUN apt-get update && apt-get -y update

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y curl

RUN apt-get install wget

RUN apt-get --assume-yes install git

RUN apt-get install -y build-essential python3.6 python3-pip python3-dev

RUN apt-get update && apt-get install -y emacs

RUN pip3 -q install pip --upgrade

RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b

RUN rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH=/miniconda/bin:${PATH}

RUN conda update -y conda

RUN conda config --add channels conda-forge

RUN mkdir src

WORKDIR src/

COPY . .

# RUN conda install -c anaconda -y python=3.7.2 https://anaconda.org/michaeljwilson/desi-high/badges/version.svg
# RUN conda install -c anaconda -y pip tensorflow-mkl=1.9.0 keras=2.2.2 

# RUN pip3 install jupyter
# RUN python3 module.py

RUN wget https://raw.githubusercontent.com/michaelJwilson/DESI-HighSchool/master/environment.yml

RUN conda env create -n desihigh -f environment.yml

RUN git clone https://github.com/michaelJwilson/DESI-HighSchool src/

WORKDIR src/

# Add Tini. Tini operates as a process subreaper for jupyter. This prevents kernel crashes.
ENV TINI_VERSION v0.6.0

ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini

RUN chmod +x /usr/bin/tini

ENTRYPOINT ["/usr/bin/tini", "--"]

CMD ["conda", "activate", "desihigh"]

CMD ["conda", "install", "jupyter"]

CMD ["git", "clone", "https://github.com/michaelJwilson/DESI-HighSchool"]

# CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]