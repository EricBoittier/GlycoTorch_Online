#FROM tiangolo/uwsgi-nginx-flask:python3.8-alpine
#
#RUN apk update
#RUN apk add make automake gcc g++ subversion python3-dev
#
#ENV STATIC_URL /static
#ENV STATIC_PATH /var/www/app/static
#COPY ./requirements.txt /var/www/requirements.txt
#COPY ./app/uwsgi.ini /app/uwsgi.ini
#
#FROM arm64v8/gcc:4.9
#FROM arm64v8/python:3

#FROM --platform=linux/amd64 arm64v8/gcc
#python:3.9-slim
#FROM  tiangolo/uwsgi-nginx-flask:python3.8-alpine
#RUN apk update
#RUN apk add make automake gcc g++ subversion python3-dev

FROM arm64v8/pypy:3
ENV STATIC_URL /static
ENV STATIC_PATH /var/www/app/static
COPY .. .
COPY requirements.txt /var/www/requirements.txt
#COPY ./uwsgi.ini /app/uwsgi.ini

RUN pip install -r /var/www/requirements.txt
