FROM ubuntu:20.04

RUN apt-get update && apt-get install -y python3 \
 python3-pip \
 python3-dev \
 python3-setuptools


RUN pip3 install --upgrade pip

RUN mkdir D2H2

COPY requirements.txt /D2H2

WORKDIR /D2H2

RUN pip3 install -r requirements.txt

COPY . .

WORKDIR /D2H2/app

EXPOSE 5000

CMD gunicorn --bind 0.0.0.0:5000 app:app
