FROM julia:1.7.0
RUN apt-get -y update
RUN apt-get -y install git

# Build
# docker build -t dpmneal3 .

# Run
# sudo sudo docker run --rm -di dpmneal3