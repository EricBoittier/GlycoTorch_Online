#FROM arm64v8/gcc:4.9
#FROM arm64v8/python:3

#FROM --platform=linux/amd64 arm64v8/gcc
#python:3.9-slim
#FROM  tiangolo/uwsgi-nginx-flask:python3.8-alpine
#RUN apk update
#RUN apk add make automake gcc g++ subversion python3-dev

#FROM arm64v8/pypy:3
FROM arm64v8/python:latest
ENV STATIC_URL /static
ENV STATIC_PATH /var/www/app/static
COPY ./requirements.txt /var/www/requirements.txt
#COPY ./uwsgi.ini /app/uwsgi.ini

COPY requirements.txt app.py ./

CMD [ "pwd" ]
CMD [ "ls" ]

RUN pip3 install -r requirements.txt

EXPOSE 8050

RUN pip install -r /var/www/requirements.txt

COPY . .

CMD [ "python3", "-m" , "flask", "run", "--host" ,"0.0.0.0", "--port", "8050"]

#CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0"]
