FROM python:3.9-bullseye
ENV DEBIAN_FRONTEND=noninteractive

RUN set -x \
  && echo "Preparing system..." \
  && apt-get -y update \
  && apt-get -y --no-install-recommends install \
    git \
    python3-dev \
    python3-pip \
    python3-setuptools \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir --upgrade pip

RUN set -x \
  && echo "Installing R..." \
  && apt-get -y update \
  && apt-get -y install r-base \
  && rm -rf /var/lib/apt/lists/* 

RUN mkdir D2H2

COPY requirements.txt /D2H2

WORKDIR /D2H2

RUN pip3 install -r requirements.txt

COPY . .

WORKDIR /D2H2/app

RUN Rscript setup.R

EXPOSE 5000

CMD gunicorn --bind 0.0.0.0:5000 app:app
